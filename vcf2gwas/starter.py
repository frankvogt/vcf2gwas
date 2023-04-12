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

import time
import os
import sys
import math
import itertools
import concurrent.futures
import multiprocessing as mp
from psutil import virtual_memory
import pandas as pd


pd.options.mode.chained_assignment = None

#################### Setting variables ####################

argvals = None

P = Parser(argvals)
out_dir = P.set_out_dir()
out_dir2 = os.path.join(out_dir, "Output")
os.makedirs(out_dir2, exist_ok=True)

timestamp = time.strftime("%a, %d %b %Y %H:%M", time.localtime())
timestamp2 = P.set_timestamp()
if timestamp2 == None:
    # for testing
    timestamp2 = time.strftime("%Y%m%d_%H%M%S")

start = time.perf_counter()

version = set_version_number()

shutil.rmtree(f'_vcf2gwas_temp_{timestamp2}', ignore_errors=True)
dir_temp = f'_vcf2gwas_temp_{timestamp2}'
qc_dir = os.path.join(dir_temp, f'QC_{timestamp2}')
make_dir(qc_dir)

file = open(os.path.join(dir_temp, "vcf2gwas_process_report.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_snpcount_total.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_snpcount_sig.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_sig_level.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_ind_count.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_ind_gemma.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_summary_paths.txt"), 'a')
file.close()
file = open(os.path.join(dir_temp, "vcf2gwas_gemma_fail.txt"), 'a')
file.close()


# get genotype file input
snp_file2 = P.set_geno()
if snp_file2.endswith((".vcf", "vcf.gz")) == False:
    msg = "VCF file is neither in .vcf or .vcf.gz format"
    raise ValueError(msg)
snp_file = os.path.split(snp_file2)[1]

# get phenotype / covariate input
pheno_files = P.set_pheno()
if pheno_files != None:
    p_list_temp = []
    pheno_dir = os.path.join(dir_temp, "Pheno")
    make_dir(pheno_dir)
    for pf in pheno_files:
        p_file = os.path.join(pheno_dir, os.path.split(pf)[1])
        shutil.copy(pf, p_file)
        p_list_temp.append(p_file)
    pheno_files = p_list_temp
pheno_files_temp = pheno_files
pheno_files_path = []
if pheno_files != None:
    plist = []
    for p in pheno_files:
        pheno_files_path.append(os.path.split(p)[0])
        plist.append(os.path.split(p)[1])
    pheno_files = plist
pheno = listtostring(pheno_files)
if pheno_files != None and len(pheno_files) > 1:
    switch = True
    pheno = pheno.replace(".csv", "")
    pheno = pheno.replace(" ", "_")
else:
    switch = False

covar = P.set_covar()
ascovariate = P.set_ascovariate()
transform_metric = P.set_transform()
transform_switch = False
if transform_metric != None:
    transform_switch = True
umap_n = P.set_U()
pca_n = P.set_P()
geno_pca_switch = False
umap_switch2 = False
pca_switch2 = False
if transform_switch == True:
    ascovariate = False
if umap_n == None and pca_n == None:
    ascovariate = False
if covar != None and ascovariate == True:
    msg = "Only one covariate file allowed, cannot use both -cf/--cfile and -asc/--ascovariate"
    raise SyntaxError(msg)
elif covar != None:
    if covar.lower() == "pca":
        covar = os.path.join(dir_temp, "vcf2gwas_geno_pca.csv")
        geno_pca_switch = True
elif ascovariate == True:
    if switch == True:
        msg = "Only one phenotype file is allowed when employing -asc/--ascovariate"
        raise SyntaxError(msg)
    elif umap_n != None and pca_n != None:
        msg = "Only one covariate file allowed, cannot use both -U/--UMAP and -P/--PCA in conjunction with -asc/--ascovariate"
        raise SyntaxError(msg)
    elif umap_n != None:
        umap_switch2 = True
        pheno_files_ = pheno_files[0].removesuffix(".csv")
        covar = os.path.join(pheno_files_path[0], f'{pheno_files_}_umap.csv')
    elif pca_n != None:
        pca_switch2 = True
        pheno_files_ = pheno_files[0].removesuffix(".csv")
        covar = os.path.join(pheno_files_path[0], f'{pheno_files_}_pca.csv')
if covar == None:
    covar_temp = covar
else:
    covar_temp = [covar]
if covar != None:
    covar_file = covar
    covar = os.path.split(covar)[1]

lm = P.set_lm()
gk = P.set_gk()
eigen = P.set_eigen()
lmm = P.set_lmm()
bslmm = P.set_bslmm()
model = set_model(lm, gk, eigen, lmm, bslmm)
if model == "-bslmm" and covar != None:
    msg = "GEMMA doesn't support adding Covariates when running BSLMM"
    raise SyntaxError(msg)
if model in ["-gk", "-eigen"]:
    pheno_files = None
    pheno = None
    covar = None
    switch = False
model2 = None
model2_temp = None
if model != None:
    model2 = change_model_dir_names(model)
    model2_temp = model.removeprefix("-")

if model == None:
    path = out_dir2
else:
    path = os.path.join(out_dir2, model2)
make_dir(path)

pc_prefix = set_pc_prefix(pheno, covar, timestamp2)
pc_prefix2 = set_pc_prefix(pheno, covar, "_")
if model in ["-gk", "-eigen"]:
    pc_prefix2 = set_pc_prefix(pheno, covar, "")

# configure logger
Log = Logger(pc_prefix, path)
Log.just_log(f"\nvcf2gwas v{version} \n")
Log.just_log("Initialising..\n")
Log.print_log(f'Start time: {timestamp}\n')
Log.print_log("Parsing arguments..")

# get gene file input
gene_file = P.set_gene_file()
species = None
if gene_file != None:
    gene_file = list(set(gene_file))
    gene_file, species = Starter.get_gene_file(gene_file, Log)
gene_file_path = []
gene_file2 = []
if gene_file != None:
    for gf in gene_file:
        gene_file_path.append(gf)
        gene_file2.append(os.path.split(gf)[1])
        gene_file = gene_file2
gene_thresh = P.set_gene_thresh()
# file checking
Log.print_log(f'Genotype file: {snp_file}')
if model in ["-lm", "-lmm", "-bslmm"] and pheno_files == None:
    msg = "No phenotype file specified"
    raise_error(SyntaxError, msg, Log)
if pheno_files != None:
    Log.print_log(f'Phenotype file(s): {listtostring(pheno_files, ", ")}')
if covar != None:
    if geno_pca_switch == True:
        Log.print_log(f'Covariates: principal components of genotype file')
    elif umap_switch2 == True:
        Log.print_log(f'Covariates: UMAP embeddings of phenotype file')
    elif pca_switch2 == True:
        Log.print_log(f'Covariates: principal components of phenotype file')
    else:
        Log.print_log(f'Covariate file: {covar}')
if gene_file != None:
    for gf, s in zip(gene_file, species):
        if s == "":
            Log.print_log(f'Gene comparison file: {gf}')
        else:
            Log.print_log(f'Gene comparison with: {s}')

pheno2 = None
covar2 = None
check_files(snp_file2, gene_file, gene_file_path)
if pheno_files != None:
    x = 0
    for pheno_file in pheno_files:
        pheno2 = os.path.join(pheno_files_path[x], pheno_file)
        check_files2(pheno_file, pheno2)
        x += 1
if covar != None:
    covar2 = covar_file
    if geno_pca_switch == False and umap_switch2 == False and pca_switch2 == False:
        check_files3(covar, covar2)

# get phenotype / covariate selection input
X = P.get_phenotypes()
if model in ["-gk", "-eigen"]:
    X = None
X_org = X
if pheno_files_temp != None and X != None:
    X = pheno_switcher(pheno_files_temp, X)
Y = P.get_covariates()
if geno_pca_switch == True:
    if Y == None:
        Y = ["2"]
        Log.print_log("Info: no number of principal components specified, now using 2 PCs for the analysis")
    n_pca = int(listtostring(Y))
    Y = None
if covar_temp == None:
    Y = None
Y_org = Y
if covar_temp != None and Y != None:
    Y = pheno_switcher(covar_temp, Y)
A = P.set_A()
B = P.set_B()
if umap_switch2 == True or pca_switch2 == True:
    B = True
if geno_pca_switch == True:
    B = True
if A == True and X != None:
    X = None
    msg = "Option 'allphenotypes' will overwrite your phenotype selection"
    raise_warning(SyntaxWarning, msg, Log)
if B == True and Y != None:
    Y = None
    msg = "Option 'allcovariates' will overwrite your covariate selection"
    raise_warning(SyntaxWarning, msg, Log)

#check for duplicate selection
x_test = []
if X != None and Y != None:
    x_test = set(X) & set(Y)
    x_test = list(x_test)
if pheno == covar and x_test != []:
    msg = f'The same data was selected as phenotype and covariate (column {listtostring(x_test, ", ")} in "{pheno}"). This will cause issues during the analysis with GEMMA'
    raise_error(ValueError, msg, Log)

# get more variables
n_top = P.set_n_top()
chr = P.set_chr()
n = set_n(lm, gk, lmm, bslmm)
filename = P.set_filename()
if filename != None:
    filename = os.path.split(filename)[1]
else:
    if model == "-eigen":
        msg = "No relatedness matrix file specified to carry out the eigen-decomposition"
        raise_error(SyntaxError, msg, Log)
min_af = P.set_q()
pca = P.set_pca()
keep = P.set_keep()
multi = P.set_multi()
if multi == True:
    switch = True
seed = P.set_seed()
sigval = P.set_sigval()
nolabel = P.set_nolabel()
noqc = P.set_noqc()
noplot = P.set_noplot()
burn = P.set_burn()
sampling = P.set_sampling()
snpmax = P.set_snpmax()

Log.print_log("\nArguments parsed successfully")

# file preparation (and get remaining variables)
Log.print_log("\nPreparing files\n")

memory = P.set_memory()
memory_org = memory
memory_total = int(((virtual_memory().total/1e9)//2)*1e3)
if memory_total < memory:
    memory = memory_total
memory2 = memory
if memory < 1000:
    msg = "Not enough memory available (at least 1000 MB required)"
    raise_error(MemoryError, msg, Log)

threads = P.set_threads()
os.environ['NUMEXPR_MAX_THREADS'] = str(threads)
threads_org = threads
cpu = mp.cpu_count()-1
if cpu == 0:
    cpu = 1
if cpu < threads:
    threads = cpu
if threads == 0:
    msg = "No logical cores available!"
    raise_error(ChildProcessError, msg, Log)

umap_switch = False
if ascovariate == False:
    if umap_n != None:
        Starter.check_vals("UMAP", umap_n, 1, 5, Log)
        umap_switch = True
umapmetric = P.set_umapmetric()

pca_switch = False
if ascovariate == False:
    if pca_n != None:
        Starter.check_vals("PCA", pca_n, 2, 10, Log)
        pca_switch = True

# check model / phenotype / genotype selection
if model == None:
    msg = "No model specified for GEMMA analysis"
    raise_error(SyntaxError, msg, Log)
if model not in ["-gk", "-eigen"]:
    if A == False:
        if umap_switch == False and pca_switch == False:
            if X == None:
                msg = "No phenotypes specified for GEMMA analysis"
                raise_error(SyntaxError, msg, Log)
    elif covar2 != None:
        if B == False:
            if geno_pca_switch == False and umap_switch2 == False and pca_switch2 == False:
                if Y == None:
                    msg = "No covariates specified for GEMMA analysis"
                    raise_error(SyntaxError, msg, Log)
    if A == False and B == False:
        if umap_switch == False and pca_switch == False:
            if geno_pca_switch == False:
                if X == None and Y == None:
                    msg = "No phenotypes and covariates specified for GEMMA analysis"
                    raise_error(SyntaxError, msg, Log)

pheno_list = []
threads_list = []
l = None
pheno_files2 = []

analysis_num = 0
part_str = "part"

if umap_switch == True:
    X = list(range(umap_n))
    X = [i+1 for i in X]
if pca_switch == True:
    X = list(range(pca_n))
    X = [i+1 for i in X]
if umap_switch == True and pca_switch == True:
    switch = True
if umap_switch == False and pca_switch == False:
    pheno_files2 = pheno_files
if X == None:
    if model not in ["-gk", "-eigen"]:
        X = []

#################### Split up files ####################

if pheno_files != None:

    if umap_switch == True and pca_switch == True:
        num = len(pheno_files)*2
    else:
        num = len(pheno_files)

    if num > threads:
        if umap_switch == True and pca_switch == True:
            msg = "Not enough logical cores available to analyze the input files when performing both UMAP and PCA"
            raise_error(ChildProcessError, msg, Log)
        else:
            msg = "Not enough logical cores available to analyze all phenotype input files"
            raise_error(ChildProcessError, msg, Log)

    if num > memory/1e3:
        if umap_switch == True and pca_switch == True:
            msg = "Not enough memory available to analyze the input files when performing both UMAP and PCA"
            raise_error(MemoryError, msg, Log)
        else:
            msg = "Not enough memory available to analyze all phenotype input files"
            raise_error(MemoryError, msg, Log)

    x = 0
    for pheno_file in pheno_files:

        pheno_path = pheno_files_path[x]
        pheno_file2 = os.path.join(pheno_path, pheno_file)
        df = Processing.load_pheno(pheno_file2)
        l = len(df.columns)

        #check if columns exist
        if umap_switch and pca_switch == False:
            if set(X).issubset([i+1 for i in list(range(l))]) == False:
                msg = "The selected phenotype data does not exist in the phenotype file"
                raise_error(ValueError, msg, Log)
        if covar2 != None:
            if Y != None:
                if geno_pca_switch == False and umap_switch2 == False and pca_switch2 == False:
                    df_covar = Processing.load_pheno(covar2)
                    l_covar = len(df_covar.columns)
                    if set(Y).issubset([i+1 for i in list(range(l_covar))]) == False:
                        msg = "The selected covariate data does not exist in the covariate file"
                        raise_error(ValueError, msg, Log)

        if transform_switch == True:
            Log.print_log(f"Performing {transform_metric} transformation on {pheno_file}..")
            pheno_file_trans = pheno_file
            pheno_file_trans = pheno_file_trans.removesuffix(".csv")
            pheno_file_trans = f'{pheno_file_trans}_transformed.csv'
            df = pheno_transformation.transform(df, pheno_file_trans, pheno_path, method=transform_metric)
            Log.print_log(f'Saved as "{pheno_file_trans}" temporarily in {pheno_path}\nTransformation successful\n')
            pheno_file = pheno_file_trans

        if umap_switch == True or umap_switch2 == True:
            if l < umap_n:
                msg = f'Not enough phenotypes in {pheno_file} to perform UMAP (at least {umap_n} required)'
                raise_error(ValueError, msg, Log)
            else:
                Log.print_log(f"Performing UMAP dimensionality reduction on {pheno_file}..")
                pheno_file3 = pheno_file
                pheno_file3 = pheno_file3.removesuffix(".csv")
                pheno_file3 = f'{pheno_file3}_umap.csv'
                if umap_switch == True:
                    pheno_files2.append(pheno_file3)
                Starter.umap_calc(df, pheno_file3, umap_n, seed, pheno_path, umapmetric, Log)
                Log.print_log(f'Saved as "{pheno_file3}" temporarily in {pheno_path}\nUMAP calculated successful\n')
        
        if pca_switch == True or pca_switch2 == True:
            if l < pca_n:
                msg = f'Not enough phenotypes in {pheno_file} to perform PCA (at least {pca_n} required)'
                raise_error(ValueError, msg, Log)
            else:
                Log.print_log(f"Performing PCA dimensionality reduction on {pheno_file}..")
                pheno_file4 = pheno_file
                pheno_file4 = pheno_file4.removesuffix(".csv")
                pheno_file4 = f'{pheno_file4}_pca.csv'
                if pca_switch == True:
                    pheno_files2.append(pheno_file4)
                Starter.pca_calc(df, pheno_file4, pca_n, pheno_path, Log)
                Log.print_log(f'Saved as "{pheno_file4}" temporarily in {pheno_path}\nPCA calculated successful\n')

        if umap_switch2 == True:
            if l < umap_n:
                msg = f'Not enough phenotypes in {pheno_file} to perform UMAP (at least {umap_n} required)'
                raise_error(ValueError, msg, Log)

        if pca_switch2 == True:
            if l < pca_n:
                msg = f'Error: not enough phenotypes in {pheno_file} to perform PCA (at least {pca_n} required)'
                raise_error(ValueError, msg, Log)

        x += 1

    x = 0
    for pheno_file in pheno_files2:
        Log.print_log(f"Checking {pheno_file}..")
        if umap_switch == True and pca_switch == True:
            pheno_path = pheno_files_path[x//2]
        else:
            pheno_path = pheno_files_path[x]
        pheno_file2 = os.path.join(pheno_path, pheno_file)
        df = Processing.load_pheno(pheno_file2)
        part_str = Starter.get_part_str(pheno_file, part_str)

        if A == False and len(X) > len(df.columns):
            msg = f'More phenotypes selected than present in {pheno_file}!\nUsing all available phenotypes..'
            raise_warning(UserWarning, msg, Log)
            X = list(range(len(df.columns)))
            X = [i+1 for i in X]
            
        if A == True:
            l = len(df.columns)
            X = list(range(len(df.columns)))
            X = [i+1 for i in X]
        elif X != []:
            l = len(X)
        else:
            l = 1

        #QC
        if noqc == False:
            try:
                QC.pheno_QC(df, X, qc_dir)
                Log.print_log("Phenotype distribution(s) successfully plotted")
            except Exception as e:
                msg = e
                raise_error(RuntimeError, msg, Log)

            
        if switch == True:
            rest = threads%len(pheno_files2)
            threads2 = threads//len(pheno_files2)
            memory2 = memory//len(pheno_files2)
            pheno_list.append(pheno_file)

        elif l == 1:
            memory2 = memory
            pheno_list.append(pheno_file)

        elif l <= threads and l <= memory/1e3:
            rest = threads%l
            threads2 = threads//l
            memory2 = memory//l
            Starter.split_phenofile1(X, df, pheno_file, pheno_list, pheno_path, part_str)
            Log.print_log("Phenotype file split up successful")

        else:
            threads_temp = threads
            if l <= threads:
                threads_temp = math.floor(memory/1e3)
            rest = 0
            threads2 = 1
            memory2 = memory//threads_temp
            col_dict = {}
            rest2 = l%threads_temp
            threads3 = l//threads_temp
            Starter.split_phenofile2(threads_temp, threads3, rest2, col_dict, X, df, pheno_file, pheno_list, pheno_path, part_str)
            Log.print_log("Phenotype file split up successful")

        x += 1
        if umap_switch == False and pca_switch == False:
            analysis_num += l

if umap_switch == True and pca_switch == True:
    analysis_num += (umap_n + pca_n)*(len(pheno_files2)//2)

if memory2 < 1000:
    msg = "Memory might not be sufficient to carry out analysis!"
    raise_warning(ResourceWarning, msg, Log)
if memory2 == 0:
    msg = "Memory not sufficient to carry out analysis!"
    raise_error(MemoryError, msg, Log)

#################### compressing, indexing and filtering VCF file ####################

if snp_file.endswith(".vcf"):
    Log.print_log("\nCompressing VCF file..")
    timer = time.perf_counter()
    snp_file2 = Converter.compress_snp_file(snp_file2)
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    snp_file = f'{snp_file}.gz'
    Log.print_log(f'VCF file successfully compressed (Duration: {runtime_format(timer_total)})')

try:
    snp_prefix = snp_file.removesuffix(".vcf.gz")
except Exception:
    snp_prefix = snp_file.removesuffix(".gz")

#index VCF
Log.print_log("\nIndexing VCF file..")
if not os.path.exists(f"{snp_file2}.csi"):
    timer = time.perf_counter()
    Converter.index_vcf(snp_file2)
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    Log.print_log(f'VCF file successfully indexed (Duration: {runtime_format(timer_total)})')
else:
    Log.print_log(f'VCF file already indexed')

chr2, chr_num = Converter.check_chrom(snp_file2, chr)
chr_list = chr2
chr3 = listtostring(chr2, ", ")
chr2 = listtostring(chr2, ",")
if chr == None:
    chr2 = None

snp_file_path = os.path.split(snp_file2)[0]
temp_dir = os.path.join(dir_temp, "VCF")
make_dir(temp_dir)
#dir_check = make_dir(temp_dir)
snp_file2_org = snp_file2
snp_file2 = os.path.join(temp_dir, snp_file)

chrom, chrom_list = Converter.set_chrom(snp_file2_org, switch=False)
if noqc == False:
    Log.print_log("\nStarting genotype Quality Control..")
    timer = time.perf_counter()
    QC.geno_QC(snp_file2_org, dir_temp, qc_dir, chrom_list, Log)
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    Log.print_log(f'Quality control successful (Duration: {runtime_format(timer_total)})')

Log.print_log("\nFiltering SNPs..")
timer = time.perf_counter()
Converter.filter_snps(min_af, snp_file2_org, snp_file2, chr2)
timer_end = time.perf_counter()
timer_total = round(timer_end - timer, 2)
Log.print_log(f'SNPs successfully filtered (Duration: {runtime_format(timer_total)})')
os.remove(f'{snp_file2_org}.csi')

if geno_pca_switch == True:
    Log.print_log("\nExtracting principal components from VCF file..")
    timer = time.perf_counter()
    list1 = Processing.process_snp_file(snp_file2)
    covar_cols = Processing.pca_analysis2(snp_file2, n_pca, memory, threads, chrom, list1, dir_temp)
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    Log.print_log(f'PCA successful (Duration: {runtime_format(timer_total)})')

#################### Prepare commands for analysis.py ####################

args_list = []
args = sys.argv[1:]
# for testing
if args == []:
    args = argvals
args = Starter.delete_string(args, ['--timestamp'])
input_str = f'vcf2gwas {listtostring(args)}'
args = Starter.delete_string(args, ['-v', '--vcf', '-T', '--threads', '-M', '--memory'])
if covar == None:
    args = Starter.delete_string(args, ['-cf', '--cfile', '-c', '--covar'])
args.insert(0, snp_file2)
args.insert(0, "--vcf")
args.insert(0, timestamp2)
args.insert(0, "--timestamp")
args.insert(0, 'python3.9')
args.insert(1, os.path.join(os.path.dirname(__file__), 'analysis.py'))
args.insert(2, '--memory')
args.insert(3, str(memory2))

if l == None:
    Starter.edit_args3(args, threads, args_list)

elif switch == True:
    Starter.adjust_threads(pheno_list, threads2, rest, threads_list)
    args = Starter.delete_string(args, ['-pf', '--pfile'])
    if umap_switch == True or pca_switch == True:
        args = Starter.delete_string(args, ['-p', '--pheno'])
    Starter.edit_args1(pheno_list, args, args_list, threads_list, umap_switch, pca_switch, A, pheno_files_path)

elif l != 1:
    Starter.adjust_threads(pheno_list, threads2, rest, threads_list)
    args = Starter.delete_string(args, ['-p', '--pheno'])
    Starter.edit_args2(pheno_list, args, args_list, threads_list, pheno, A, pheno_files_path)
    if umap_switch == True:
        Log.print_log(f'\nInfo:\nAfter reducing dimensions of {pheno} via UMAP, it has been split up in {len(pheno_list)} parts in order to ensure maximum efficiency')
    else:
        Log.print_log(f'\nInfo:\n{pheno} has been split up in {len(pheno_list)} parts in order to ensure maximum efficiency')
else:
    Starter.edit_args3(args, threads, args_list)

if covar != None:
    args_list2 = []
    for args in args_list:
        args = Starter.delete_string(args, ['-cf', '--cfile'])
        if B == True:
            args = Starter.delete_string(args, ['-c', '--covar'])
            args.insert(4, "--allcovariates")
        args.insert(4, covar2)
        args.insert(4, "--cfile")
        args_list2.append(args)
    args_list = args_list2

if model in ["-gk", "-eigen"]:
    args_list2 = []
    for args in args_list:
        args = Starter.delete_string(args, ['-pf', '--pfile', '-cf', '--cfile', '-p', '--pheno', '-c', '--covar'])
        args_list2.append(args)
    args_list = args_list2

log_path2 = os.path.join(path, "Logs", f'logs{pc_prefix2}_{snp_prefix}_{timestamp2}')
make_dir(log_path2)

Log.print_log("\nFile preparations completed")

#################### Run analysis.py for analysis ####################

Log.print_log("\nStarting analysis..")
with concurrent.futures.ProcessPoolExecutor(mp_context=mp.get_context('fork'), max_workers=threads_org) as executor:
        executor.map(Starter.run_vcf2gwas, args_list, itertools.repeat(dir_temp))
Starter.check_return_codes(Log, dir_temp)

#################### summary and clean-up ####################

snp_file2 = snp_file2_org

# summarizer and gene comparison
if model not in ("-gk", "-eigen", None):
    path2_temp = os.path.join(path, "Summary", "temp")
    make_dir(path2_temp)
    path2 = os.path.join(path, "Summary", f'summary{pc_prefix2}_{snp_prefix}_{timestamp2}')
    os.rename(path2_temp, path2)
    path3 = os.path.join(path2, "top_SNPs")
    prefix_list = []
    for pheno in pheno_list:
        pc_prefix3 = set_pc_prefix(pheno, covar, "_")
        prefix_list.append(pc_prefix3)
    filenames, str_list = Summary.summarizer(path3, path2, pc_prefix3, snp_prefix, n_top, Log, prefix_list)
    temp, file_dict = Summary.ind_summary(path2, filenames, str_list, dir_temp)
    filenames = temp[0]
    str_list = temp[1]
    filenames2 = []
    name_list = []
    name_list2 = []
    if gene_file != None:
        filenames2, temp = Summary.gene_compare(filenames, str_list, gene_file, gene_file_path, gene_thresh, path2, snp_prefix, chr_list, Log)
        name_list = temp[0]
        name_list2 = temp[1]
        Summary.gene_occurrence(filenames2)
        Summary.pheno_compare_split(filenames2, file_dict, name_list, name_list2)

#move QC files
if noqc == False:
    if pheno_files != None:
        qc_path = os.path.join(path, "QC")
        make_dir(qc_path)
        shutil.move(qc_dir, qc_path)
        Log.print_log(f'\nQuality control files moved to {qc_path}')

# move split up files (and reduce "files" folder)
path5 = os.path.join(path, "Files")
pheno_temp = [f'files_{x.removesuffix(".csv")}_{snp_prefix}' for x in pheno_list]
if switch == False:
    if len(pheno_list) > 1:
        for file in pheno_list:
            os.remove(os.path.join(pheno_files_path[0], file))
        c = 1
        for pheno_t in pheno_temp:                
            for folder in os.listdir(path5):
                if pheno_t == folder:
                    folder2 = folder.replace(f".{part_str}{c}", "")
                    folder2 = f'{folder2}_{timestamp2}'
                    shutil.move(os.path.join(path5, folder), os.path.join(path5, folder2))
                    path6 = os.path.join(path5, folder2)
                    for file in os.listdir(path6):
                        for string in [".log.txt", ".cXX.txt", ".log", ".bim", ".bed", ".nosex", ".covariates.txt"]:
                            if file.endswith(string):
                                os.rename(os.path.join(path6, file), os.path.join(path6, file.replace(f".{part_str}{c}", "")))
                    for folder3 in os.listdir(path5):
                        for pheno in pheno_temp:
                            if pheno == folder3:
                                for file in os.listdir(os.path.join(path5, folder3)):
                                    if file.endswith(".fam"):
                                        shutil.move(os.path.join(path5, folder3, file), os.path.join(path5, folder2, file))
                                shutil.rmtree(os.path.join(path5, folder3))
            c += 1
    else:
        for folder in os.listdir(path5):
            if model in ["-gk", "-eigen"]:
                if f'files_{snp_prefix}' == folder:
                    shutil.move(os.path.join(path5, folder), os.path.join(path5, f'{folder}_{timestamp2}'))
            elif pheno_temp[0] == folder:
                os.rename(os.path.join(path5, folder), os.path.join(path5, f'{folder}_{timestamp2}'))
elif switch == True:
    for folder in os.listdir(path5):
        for pheno in pheno_temp:
            if pheno == folder:
                os.rename(os.path.join(path5, folder), os.path.join(path5, f'{folder}_{timestamp2}'))
if geno_pca_switch == True:
    for folder in os.listdir(path5):
        if folder.endswith(timestamp2):
            shutil.copy(covar2, os.path.join(path5, folder, f'{snp_prefix}_PCA.csv'))

# move umap/pca files
if umap_switch == True or pca_switch ==True or umap_switch2 == True or pca_switch2 == True:
    switch_names = ["UMAP", "PCA"]
    x = 0
    for switch in [umap_switch, pca_switch, umap_switch2, pca_switch2]:
        if switch == True:
            path4 = os.path.join(path, switch_names[x%2], f'{switch_names[x%2]}_{timestamp2}')
            make_dir(path4)
            y = 0 
            for pheno_file in pheno_files:
                pheno_path = pheno_files_path[y]
                if umap_switch == True and pca_switch == True:
                    pheno_path = pheno_files_path[y//2]
                for files in os.listdir(pheno_path):
                    if transform_switch == True:
                        if files.startswith(f'{pheno_file.removesuffix(".csv")}_transformed_{switch_names[x%2].lower()}'):
                            shutil.move(os.path.join(pheno_path, files), os.path.join(path4, files))
                    else:
                        if files.startswith(f'{pheno_file.removesuffix(".csv")}_{switch_names[x%2].lower()}'):
                            shutil.move(os.path.join(pheno_path, files), os.path.join(path4, files))
                y += 1
            Log.print_log(f'\nMoved {switch_names[x%2]} files to {path4}')
        x += 1

# move transformed pheno files
if transform_switch == True:
    path5 = os.path.join(path, "Transformation", f'{"Transformation"}_{timestamp2}')
    make_dir(path5)
    y = 0 
    for pheno_file in pheno_files:
        pheno_path = pheno_files_path[y]
        for files in os.listdir(pheno_path):
            if files.startswith(f'{pheno_file.removesuffix(".csv")}_transformed'):
                shutil.move(os.path.join(pheno_path, files), os.path.join(path5, files))
        y += 1

# remove temp covariate file
Converter.remove_covar(covar)

finish = time.perf_counter()
time_total = round(finish-start, 2)
time_total = runtime_format(time_total)

Log.print_log(f'\nClean up successful \n\nvcf2gwas has been successfully completed! \nRuntime: {time_total}\n')

#get variables for summary
if pheno_files != None:
    pheno_files = listtostring(pheno_files, ', ')
X_names = Y_names = ""

if X != None:
    X1 = []
    if umap_switch == True:
        X1.append(f'{umap_n} UMAP embedding(s)')
    if pca_switch == True:
        X1.append(f'{pca_n} principal components')
    if umap_switch == False and pca_switch == False:
        X, X_names = pheno_switcher2(pheno_files_temp, X_org, X)
        Summary.pheno_summary(filenames, filenames2, str_list, path2, X_names, name_list, snp_prefix)
        X_names = listtostring(X_names, ", ")
        X1 = X
    X = listtostring(X1, ', ')
    for f in os.listdir(path2):
        if "temp_summary" in f:
            os.remove(os.path.join(path2, f))
if Y != None:
    Y, Y_names = pheno_switcher2(covar_temp, Y_org, Y)
    Y_names = listtostring(Y_names, ", ")
    Y = listtostring(Y, ', ')
if geno_pca_switch == True:
    Y = len(covar_cols)
if umap_switch2 == True:
    Y = umap_n
if pca_switch2 == True:
    Y = pca_n

# print summary and move log files
snp_total, snp_sig = Starter.get_snpcounts(dir_temp)
if noplot == True:
    snp_sig = "-"
sig_level = Starter.get_count("vcf2gwas_sig_level.txt", dir_temp)
ind_count = Starter.get_count("vcf2gwas_ind_count.txt", dir_temp)
gemma_count = Starter.get_count("vcf2gwas_ind_gemma.txt", dir_temp)
failed_count, failed_list = Starter.get_gemma_fail(dir_temp)
Log.summary(
    snp_file, pheno_files, covar, X, Y, model2_temp, n, filename, min_af, A, B, 
    pca, keep, memory, threads, n_top, gene_file, species, gene_thresh, multi, umap_n, umapmetric, pca_n, 
    out_dir2, analysis_num, sigval, nolabel, chr, chr3, chr_num, X_names, snp_total, snp_sig, sig_level, geno_pca_switch, burn, sampling, snpmax, noqc,
    input_str, noplot, ind_count, gemma_count, umap_switch2, pca_switch2, ascovariate, transform_metric, failed_count, failed_list
)

if model != None:
    shutil.move(os.path.join(path, f'vcf2gwas{pc_prefix}.log.txt'), os.path.join(path, f'vcf2gwas_{snp_prefix}{pc_prefix2}_{timestamp2}.log.txt'))

shutil.rmtree(dir_temp, ignore_errors=True)
