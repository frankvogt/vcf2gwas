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

from vcf2gwas.parsing import *
from vcf2gwas.utils import *
#from parsing import *
#from utils import *

import sys
import time
import concurrent.futures
import multiprocessing as mp
import itertools
import os
import shutil
import pandas as pd

#os.chdir(os.path.dirname(os.path.abspath(__file__)))

argvals = None

############################## Initialising Program ##############################

start = time.perf_counter()
# set argument parser
P = Parser(argvals)
out_dir = P.set_out_dir()
out_dir2 = os.path.join(out_dir, "Output")
timestamp2 = P.set_timestamp()
dir_temp = f'_vcf2gwas_temp_{timestamp2}'
memory = P.set_memory()
threads = P.set_threads()
os.environ['NUMEXPR_MAX_THREADS'] = str(threads)

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
    model2 = change_model_dir_names(model)

pc_prefix = set_pc_prefix(pheno_file, covar_file, "_")
if model in ["-gk", "-eigen"]:
    pc_prefix = set_pc_prefix(pheno_file, covar_file, "")

if model == None:
    path = out_dir2
else:
    path = os.path.join(out_dir2, model2)

# configure logger
log_path = get_log_path(path, timestamp2)
Log = Logger(pc_prefix, log_path)

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

n_top = P.set_n_top()
gene_thresh = P.set_gene_thresh()
multi = P.set_multi()
sigval = P.set_sigval()
nolabel = P.set_nolabel()
noplot = P.set_noplot()
burn = str(P.set_burn())
sampling = str(P.set_sampling())
snpmax = str(P.set_snpmax())

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
chrom, chrom_list = Converter.set_chrom(snp_file2)
subset = f'sub{pc_prefix}_{snp_prefix}'
subset_temp_1 = subset
folder = f'files{pc_prefix}_{snp_prefix}'

subset2 = os.path.join(dir_temp, f'mod_{subset}')
subset = os.path.join(dir_temp, subset)
subset_temp = subset

# make list from genotype file
Log.print_log("Checking individuals in VCF file..")
list1 = Processing.process_snp_file(snp_file2)
list1_org = list1

# make list from phenotype and/or covariate file and remove overlap in all files
diff1a = []
diff1b = []

if pheno_file == None:
    Log.print_log("No phenotype file specified")
    pheno_subset1 = pd.DataFrame()
    subset_temp = snp_file2
else:
    Log.print_log("Checking individuals in phenotype file..")
    if covar_file != None:
        subset_temp_1 = f'temp_{subset_temp_1}'
        subset_temp = os.path.join(dir_temp, subset_temp_1)
    pheno = Processing.load_pheno(pheno_file2)
    list2 = Processing.pheno_index(pheno)
    diff1a = Processing.make_diff(list1, list2)
    diff2 = Processing.make_diff(list2, list1)
    pheno_subset1 = Processing.make_uniform(list1, list2, diff1a, diff2, pheno, subset_temp, snp_file2, pheno_file, "phenotype", Log)
    length1 = len(pheno_subset1.columns)
    if covar_file != None:
        subset_temp_1 = f'{subset_temp_1}.vcf.gz'
        subset_temp = os.path.join(dir_temp, subset_temp_1)
    Log.print_log(f'Removed {len(diff1a)} out of {len(list1)} genotype individuals, {len(list1)-len(diff1a)} remaining \nRemoved {len(diff2)} out of {len(list2)} phenotype individuals, {len(list2)-len(diff2)} remaining')
    list1 = [i for i in list1 if i not in diff1a]

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
        Log.print_log("Not all individuals in covariate and genotype file match")
    else:
        pheno_subset2 = Processing.make_uniform(list1, list3, diff1b, diff3, covar, subset, subset_temp, covar_file, "covariate", Log)
    length2 = len(pheno_subset2.columns)
    Log.print_log(f'Removed {len(diff1b)} out of {len(list1)} genotype individuals, {len(list1)-len(diff1b)} remaining \nRemoved {len(diff3)} out of {len(list3)} covariate individuals, {len(list3)-len(diff3)} remaining')

if model in ("-gk", "-eigen"):
    if pheno_file == None and covar_file == None:
        shutil.copy(snp_file2, f'{subset}.vcf.gz')

diff1c = list(set(diff1a) & set(diff1b))
diff_num = len(diff1a)+len(diff1b)-len(diff1c)

if len(list1_org)-diff_num == 0:
    if keep == False:
        Converter.remove_files(subset, pheno_file, subset2, snp_file2)
    msg = "No individuals left! Check if IDs in VCF and phenotype file are of the same format"
    raise_error(ValueError, msg, Log)

file = open(os.path.join(dir_temp, "vcf2gwas_ind_count.txt"), 'a')
file.write(f'{len(list1_org)-diff_num}\n')
file.close()

Log.print_log(f'In total, removed {diff_num} out of {len(list1_org)} genotype individuals, {len(list1_org)-diff_num} remaining')
Log.print_log("Files successfully adjusted\n")

############################## Filter, make ped and bed ##############################

Log.print_log("Filtering and converting files\n")
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
        cols1 = Processing.edit_fam(fam, pheno_subset1, subset2, X, "p", "Phenotype", Log)
        cols1 = [str(i) for i in cols1]
    else:
        cols1 = []
        if X == None:
            X = []
            Log.print_log("No phenotypes were specified, continuing without")
        else:
            Processing.edit_fam(fam, pheno_subset1, subset2, X, "p", "Phenotype", Log)

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
        covar_file_name = Processing.make_covarfile(fam, pheno_subset2, subset2, Y, Log)
        Y = []
        cols2 = []
    elif B == True:
        Y = list(range(length2))
        Y = [i+1 for i in Y]
        Log.print_log("All covariates chosen")
        cols2 = Processing.edit_fam(fam, pheno_subset2, subset2, Y, "c", "Covariate", Log)
        cols2 = [str(i) for i in cols2]
    else:
        cols2 = []
        if Y == None:
            Y = []
            Log.print_log("No covariates were specified, continuing without")
        else:
            Processing.edit_fam(fam, pheno_subset2, subset2, Y, "c", "Covariate", Log)

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
top_sig = []
top_all = []

pd.set_option('use_inf_as_na', True)
pd.options.mode.chained_assignment = None

############################## GEMMA ##############################

prefix_list = []
prefix2_list = []
i_list = []
i_list2 = []
N_list = []
path_list = []

file = open(os.path.join(dir_temp, f"vcf2gwas_process_report_{pc_prefix}.txt"), 'a')
file.close()

if model == None:
    Log.print_log("GEMMA can't be executed since no model was specified!\n")
else:
    Log.print_log("Running GEMMA\n")
    Log.print_log(f'Phenotypes to analyze: {listtostring(columns, ", ")}\n')
    timer = time.perf_counter()
    # set lists of variables and make output directories
    for i in columns:

        prefix = subset2
        prefix2 = f'{i}_{os.path.split(prefix)[1]}'

        if i == "":
            prefix2 = prefix
            i = "results"

        if X == [] and Y == []:
            N = [1]
            #N = "1"
        else:
            N = [x]
            #N = str(x)
            x = x + 1
            if multi == True:
                N = concat_lists(X, Y)
                #N = listtostring(N)

        path_temp = None
        if model not in ("-gk", "-eigen"):
            path_temp = os.path.join(path, i, f'{i}_{timestamp2}')
            make_dir(path_temp)
            file = open(os.path.join(path_temp, f'Summary_{prefix2}.txt'), 'a')
            file.write(f'Individuals cleared for analyis: {len(list1_org)-diff_num}\n')
            file.close()
        
        prefix_list.append(prefix)
        prefix2_list.append(prefix2)
        i_list.append(i)
        N_list.append(N)
        path_list.append(path_temp)

        if model in ("-eigen", "-lmm"):
            if filename == None:
                filename, code = Gemma.rel_matrix(prefix, Log, covar_file_name, pc_prefix)
                if code != "0":
                    Gemma.write_returncodes(code, pc_prefix, dir_temp)

    # run GEMMA in parallel
    with concurrent.futures.ProcessPoolExecutor(mp_context=mp.get_context('fork'), max_workers=threads) as executor:
        executor.map(
            Gemma.run_gemma, prefix_list, prefix2_list, itertools.repeat(model), itertools.repeat(n), N_list, path_list, itertools.repeat(Log), itertools.repeat(filename), 
            itertools.repeat(filename2), itertools.repeat(pca), itertools.repeat(covar_file_name), i_list, itertools.repeat(i_list2),
            itertools.repeat(burn), itertools.repeat(sampling), itertools.repeat(snpmax), itertools.repeat(pc_prefix), itertools.repeat(dir_temp)
        )
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    Post_analysis.check_return_codes(pc_prefix, dir_temp)
    Post_analysis.get_gemma_success(i_list, prefix2_list, path_list, columns, i_list2, dir_temp)
    Log.print_log(f'\nGEMMA completed successfully (Duration: {runtime_format(timer_total)})\n')

    ############################## Processing and plotting ##############################

    Log.print_log("Analyzing GEMMA results\n")
    for (top_ten, top_sig, top_all, Log, model, n, prefix2, path_temp, n_top, i, sigval, nolabel, noplot) in zip(
        itertools.repeat(top_ten), itertools.repeat(top_sig), itertools.repeat(top_all), itertools.repeat(Log), itertools.repeat(model), itertools.repeat(n), prefix2_list, path_list, 
        itertools.repeat(n_top), i_list, itertools.repeat(sigval), itertools.repeat(nolabel), itertools.repeat(noplot)
        ):
        Post_analysis.run_postprocessing(top_ten, top_sig, top_all, Log, model, n, prefix2, path_temp, n_top, i, sigval, nolabel, noplot, dir_temp)
    Log.print_log("Analysis of GEMMA results completed successfully\n")

############################## Summary and Clean up ##############################

Log.print_log("Starting clean-up\n")

if model not in ("-gk", "-eigen", None):
    path2 = os.path.join(path, "Summary", "temp","top_SNPs")
    make_dir(path2)
    Post_analysis.print_top_list(top_ten, top_sig, top_all, columns, path2, pc_prefix, snp_prefix)
    Log.print_log("Top SNPs saved")

make_dir(os.path.join(path, "Files"))
path3 = os.path.join(path, "Files", folder)
make_dir(path3)

Log.print_log("Moving files..")
for files in os.listdir(dir_temp):
    if files.startswith((os.path.split(subset)[1], os.path.split(subset2)[1])):
        shutil.move(os.path.join(dir_temp, files), os.path.join(path3, files))
if model == "-gk":
    path2 = os.path.join(path, "Relatedness Matrix", f'rel_matrix_{timestamp2}')
    make_dir(path2)
    for files in os.listdir(path3):
        if files.endswith("XX.txt"):
            shutil.move(os.path.join(path3, files), os.path.join(path2, files))
if model == "-eigen":
    path2 = os.path.join(path, "Eigen Vectors", f'rel_matrix_{timestamp2}')
    make_dir(path2)
    for files in os.listdir(path3):
        if files.endswith(("eigenD.txt", "eigenU.txt")):
            shutil.move(os.path.join(path3, files), os.path.join(path2, files))

finish = time.perf_counter()
time_total = round(finish-start, 2)
time_total = runtime_format(time_total)

Log.print_log(f'Clean up successful \n\nAnalysis of {pheno_file} finished successfully\nRuntime: {time_total}\n')

move_log(model, pc_prefix, snp_prefix, timestamp2, log_path)
