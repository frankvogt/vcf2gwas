#!/usr/bin/env python

"""
Copyright (C) 2021, Frank Vogt

This file is part of vcf2gwas.

vcf2gwas is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vcf2gwas is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vcf2gwas.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
from vcf2gwas.utils import make_dir, runtime_format
from parsing import *
from utils import *

import time
import concurrent.futures
import multiprocessing as mp
import itertools
import os
import shutil
import signal

import pandas as pd
import seaborn as sns

#os.chdir(os.path.dirname(os.path.abspath(__file__)))

argvals = None

############################## Initialising Program ##############################

timestamp = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
timestamp2 = read_timestamp()
start = time.perf_counter()
# set argument parser
P = Parser(argvals)
out_dir = P.set_out_dir()
out_dir2 = os.path.join(out_dir, "output")

# set variables
pheno_file = P.set_pheno()
pheno_file_temp = pheno_file
pheno_file_path = []
if pheno_file != None:
    plist = []
    for p in pheno_file:
        pheno_file_path.append(os.path.split(p)[0])
        pheno_file_path = listtostring(pheno_file_path)
        plist.append(os.path.split(p)[1])
    pheno_file = plist
pheno_file = listtostring(pheno_file)

covar_file = P.set_covar()
covar_temp = [covar_file]
if covar_file != None:
    covar_file_path = covar_file
    covar_file = os.path.split(covar_file)[1]

lm = P.set_lm()
gk = P.set_gk()
eigen = P.set_eigen()
lmm = P.set_lmm()
bslmm = P.set_bslmm()
model = set_model(lm, gk, eigen, lmm, bslmm)
model2 = None
if model != None:
    model2 = model.removeprefix("-")

pc_prefix = set_pc_prefix(pheno_file, covar_file, "_")
if model in ["-gk", "-eigen"]:
    pc_prefix = set_pc_prefix(pheno_file, covar_file, "")

# configure logger
Log = Logger(pc_prefix, out_dir2)

Log.print_log(f'\nBeginning with analysis of {pheno_file}\n')

# set variables
snp_file2 = P.set_geno()
snp_file = os.path.split(snp_file2)[1]
X = P.get_phenotypes()
if pheno_file_temp != None and X != None:
    X = pheno_switcher(pheno_file_temp, X)
Y = P.get_covariates()
if covar_temp != None and Y != None:
    Y = pheno_switcher(covar_temp, Y)
min_af = P.set_q()
A = P.set_A()
B = P.set_B()
pca = P.set_pca()
keep = P.set_keep()
memory = P.set_memory()
threads = P.set_threads()

n_top = P.set_n_top()
#gene_file = P.set_gene_file()
gene_thresh = P.set_gene_thresh()
multi = P.set_multi()
sigval = P.set_sigval()
nolabel = P.set_nolabel()

#snp_file2 = f'input/{snp_file}'
pheno_file2 = None
covar_file2 = None

if pheno_file != None:
    pheno_file2 = os.path.join(pheno_file_path, pheno_file)
if covar_file != None:
    covar_file2 = covar_file_path

covar_file_name = None

############################## Adjust Files ##############################

Log.print_log("Preparing files\n")
Log.print_log("Checking and adjusting files..")

snp_prefix = snp_file.removesuffix(".vcf.gz")

chrom = Converter.set_chrom(snp_file2)

subset = f'sub{pc_prefix}_{snp_prefix}'
folder = f'files{pc_prefix}_{snp_prefix}'

subset2 = f'mod_{subset}'

# make list from genotype file
Log.print_log("Checking individuals in VCF file..")
list1 = Processing.process_snp_file(snp_file2)

# make list from phenotype and/or covariate file and remove overlap in all files
diff1a = []
diff1b = []

if pheno_file == None:
    Log.print_log("No phenotype file specified")
    pheno_subset1 = pd.DataFrame()
else:
    Log.print_log("Checking individuals in phenotype file..")
    pheno = Processing.load_pheno(pheno_file2)
    list2 = Processing.pheno_index(pheno)
    diff1a = Processing.make_diff(list1, list2)
    diff2 = Processing.make_diff(list2, list1)
    pheno_subset1 = Processing.make_uniform(list1, list2, diff1a, diff2, pheno, subset, snp_file2, pheno_file, "phenotype", Log)
    length1 = len(pheno_subset1.columns)
    Log.print_log(f'Removed {len(diff2)} out of {len(list2)} individuals, {(len(list2)-len(diff2))} remaining')

if covar_file == None:
    Log.print_log("No covariate file specified")
    pheno_subset2 = pd.DataFrame()
else:
    Log.print_log("Checking individuals in covariate file..")
    covar = Processing.load_pheno(covar_file2)
    list3 = Processing.pheno_index(covar)
    diff1b = Processing.make_diff(list1, list3)
    diff3 = Processing.make_diff(list3, list1)
    if set(diff1a) == set(diff1b) and len(set(diff1a)) != 0:       
        pheno_subset2 = Processing.rm_pheno(covar, diff3, covar_file)
        Log.print_log(str("Not all individuals in covariate and genotype file match"))
    else:
        pheno_subset2 = Processing.make_uniform(list1, list3, diff1b, diff3, covar, subset, snp_file2, covar_file, "covariate", Log)
    length2 = len(pheno_subset2.columns)
    Log.print_log(f'Removed {len(diff3)} out of {len(list3)} individuals, {(len(list3)-len(diff3))} remaining')

if model in ("-gk", "-eigen"):
    if pheno_file == None and covar_file == None:
        shutil.copy(snp_file2, f'{subset}.vcf.gz')

Log.print_log(f'In total, removed {(len(diff1a)+len(diff1b))} out of {len(list1)} individuals, {(len(list1)-(len(diff1a)+len(diff1b)))} remaining')
Log.print_log("Files successfully adjusted\n")

if (len(list1)-(len(diff1a)+len(diff1b))) == 0:
    if keep == False:
        Converter.remove_files(subset, pheno_file, subset2, snp_file2)
    Log.print_log("Error: No individuals left! Check if IDs in VCF and phenotype file are of the same format")
    sys.exit(1)

############################## Filter, make ped and bed ##############################

Log.print_log("Filtering and converting files\n")

#Log.print_log("Filtering SNPs..")
#Converter.filter_snps(min_af, subset, subset2)
#Log.print_log("SNPs successfully filtered")
os.rename(f'{subset}.vcf.gz', f'{subset2}.vcf.gz')

Log.print_log("Converting to PLINK BED..")
timer = time.perf_counter()
Converter.make_bed(subset2, chrom, memory, threads, list1)
timer_end = time.perf_counter()
timer_total = round(timer_end - timer, 2)
Log.print_log(f'Successfully converted to PLINK BED (Duration: {runtime_format(timer_total)})\n')

#remove old files
if keep == False:
    Converter.remove_files(subset, pheno_file, subset2, snp_file2)

if pca != None:
    Log.print_log("Performing principal component analysis..")
    Processing.pca_analysis(subset2, pca, memory, threads, chrom)
    Log.print_log("Principal component analysis successfully completed\n")

############################## add phenotypes/covariates to .fam file ##############################

Log.print_log("Adding phenotypes/covariates to .fam file\n")

Log.print_log("Editing .fam file..")
fam = Processing.prepare_fam(subset2)

if pheno_file == None:
    X = []
    cols1 = []
    Log.print_log("No phenotype file specified, continuing without")
else:
    if A == True:
        X = list(range(length1))
        X = [i+1 for i in X]
        Log.print_log("All phenotypes chosen")
        cols1 = Processing.edit_fam(fam, pheno_subset1, subset2, X, "p", "Phenotype", Log, model, model2, pc_prefix)
        cols1 = [str(i) for i in cols1]
    else:
        cols1 = []
        if X == None:
            X = []
            Log.print_log("No phenotypes were specified, continuing without")
        else:
            Processing.edit_fam(fam, pheno_subset1, subset2, X, "p", "Phenotype", Log, model, model2, pc_prefix)

if covar_file == None:
    Y = []
    cols2 = []
    if model == "-lmm":
        Log.print_log("No covariate file specified, continuing without")
else:
    if model in ["-lm", "-lmm"]:
        if B == True:
            Y = list(range(length2))
            Y = [i+1 for i in Y]
        covar_file_name = Processing.make_covarfile(fam, pheno_subset2, subset2, Y)
        Y = []
        cols2 = []
    elif B == True:
        Y = list(range(length2))
        Y = [i+1 for i in Y]
        Log.print_log("All covariates chosen")
        cols2 = Processing.edit_fam(fam, pheno_subset2, subset2, Y, "c", "Covariate", Log, model, model2, pc_prefix)
        cols2 = [str(i) for i in cols2]
    else:
        cols2 = []
        if Y == None:
            Y = []
            Log.print_log("No covariates were specified, continuing without")
        else:
            Processing.edit_fam(fam, pheno_subset2, subset2, Y, "c", "Covariate", Log, model, model2, pc_prefix)

if model in ("-gk", "-eigen"):
    if pheno_file == None and covar_file == None:
        fam['results'] = 1
        fam.to_csv(f'{subset2}.fam', sep=' ', header=False)

Log.print_log("Editing .fam file successful\n")

############################## Initialising GEMMA ##############################

Log.print_log("Initialising GEMMA\n")
# set variables
if cols1 == []:
        cols1 = set_cols(cols1, X, pheno_subset1)
if cols2 == []:
        cols2 = set_cols(cols2, Y, pheno_subset2)

if multi == True:
    columns = [listtostring(cols1 + cols2, '+')]
else:
    columns = cols1 + cols2
if model == "-gk" or model == "-eigen":
    columns = ['']
if X == [] and Y == []:
    columns = ['']

n = set_n(lm, gk, lmm, bslmm)
filename = P.set_filename()
#if filename != None:
#    filename = f'input/{filename}'
filename2 = None
if pca != None:
    filename = f'{subset2}.eigenvec'
    filename2 = f'{subset2}.eigenval'

x = 1
top_ten = []
sns.set_theme(style="white")
pd.set_option('use_inf_as_na', True)
pd.options.mode.chained_assignment = None

############################## GEMMA ##############################

prefix_list = []
prefix2_list = []
i_list = []
N_list = []
path_list = []

if model == None:
    Log.print_log("GEMMA can't be executed since no model was specified!\n")
else:
    Log.print_log("Running GEMMA\n")
    Log.print_log(f'Phenotypes to analyze: {listtostring(columns, ", ")}\n')
    timer = time.perf_counter()
    # set lists of variables and make output directories
    for i in columns:

        prefix = subset2
        prefix2 = f'{i}_{prefix}'
        if i == "":
            prefix2 = prefix
            i = "results"

        if X == [] and Y == []:
            N = "1"
        else:
            N = str(x)
            x = x + 1
            if multi == True:
                N = concat_lists(X, Y)
                N = listtostring(N)

        make_dir(os.path.join(out_dir2, model2))
        path = None
        if model not in ("-gk", "-eigen"):
            path = os.path.join(out_dir2, model2, i, f'{i}_{timestamp2}')
            make_dir(path)
        
        prefix_list.append(prefix)
        prefix2_list.append(prefix2)
        i_list.append(i)
        N_list.append(N)
        path_list.append(path)

        if model in ("-eigen", "-lmm"):
            if filename == None:
                filename = Gemma.rel_matrix(prefix, Log)

    # run GEMMA in parallel
    with concurrent.futures.ProcessPoolExecutor(mp_context=mp.get_context('fork'), max_workers=threads) as executor:
        executor.map(Gemma.run_gemma, prefix_list, prefix2_list, itertools.repeat(model), itertools.repeat(n), N_list, path_list, itertools.repeat(Log), itertools.repeat(filename), itertools.repeat(filename2), itertools.repeat(pca), itertools.repeat(covar_file_name), i_list)
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    Post_analysis.check_return_codes()
    Log.print_log(f'\nGEMMA completed successfully (Duration: {runtime_format(timer_total)})\n')

    ############################## Processing and plotting ##############################

    Log.print_log("Analyzing GEMMA results\n")
    for (top_ten, Log, model, n, prefix2, path, n_top, i, sigval, nolabel) in zip(itertools.repeat(top_ten), itertools.repeat(Log), itertools.repeat(model), itertools.repeat(n), prefix2_list, path_list, itertools.repeat(n_top), i_list, itertools.repeat(sigval), itertools.repeat(nolabel)):
        Post_analysis.run_postprocessing(top_ten, Log, model, n, prefix2, path, n_top, i, sigval, nolabel)
    Log.print_log("Analysis of GEMMA results completed successfully\n")

############################## Summary and Clean up ##############################

Log.print_log("Starting clean-up\n")

if model == None:
    path = out_dir2
else:
    path = os.path.join(out_dir2, model2)

if model not in ("-gk", "-eigen", None):
    #make_dir(os.path.join(path, "summary"))
    path2 = os.path.join(path, "summary", "temp","top_SNPs")
    make_dir(path2)
    Post_analysis.print_top_list(top_ten, columns, path2, pc_prefix, snp_prefix)
    Log.print_log("Top SNPs saved")

make_dir(os.path.join(path, "files"))
path3 = os.path.join(path, "files", folder)
make_dir(path3)

Log.print_log("Moving files..")
for files in os.listdir():
    if files.startswith((f'sub{pc_prefix}', subset2)):
        shutil.move(files, os.path.join(path3, files))
if model == "-gk":
    path2 = os.path.join(path, "rel_matrix", f'rel_matrix_{timestamp2}')
    make_dir(path2)
    for files in os.listdir(path3):
        if files.endswith("XX.txt"):
            shutil.move(os.path.join(path3, files), os.path.join(path2, files))
if model == "-eigen":
    path2 = os.path.join(path, "eigen_v", f'rel_matrix_{timestamp2}')
    make_dir(path2)
    for files in os.listdir(path3):
        if files.endswith(("eigenD.txt", "eigenU.txt")):
            shutil.move(os.path.join(path3, files), os.path.join(path2, files))

finish = time.perf_counter()
time_total = round(finish-start, 2)
time_total = runtime_format(time_total)

Log.print_log(f'Clean up successful \n\nAnalysis of {pheno_file} finished successfully\nRuntime: {time_total}\n')

move_log(model, model2, pc_prefix, snp_prefix, out_dir2, timestamp2)

sys.exit(0)