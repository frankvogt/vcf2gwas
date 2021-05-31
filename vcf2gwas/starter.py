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

from vcf2gwas.utils import Starter
from parsing import *
from utils import *

import subprocess
import time
import os
import sys
import concurrent.futures
import multiprocessing as mp

from psutil import virtual_memory
import pandas as pd


#os.chdir(os.path.dirname(os.path.abspath(__file__)))
pd.options.mode.chained_assignment = None

#################### Setting variables ####################

argvals = None

P = Parser(argvals)
out_dir = P.set_out_dir()
out_dir2 = os.path.join(out_dir, "output")
os.makedirs(out_dir2, exist_ok=True)

timestamp = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
start = time.perf_counter()

snp_file2 = P.set_geno()
snp_file = os.path.split(snp_file2)[1]

pheno_files = P.set_pheno()
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

#covar_files = P.set_covar()
#covar = listtostring(covar_files)
covar = P.set_covar()
if covar != None:
    covar_file = covar
    covar = os.path.split(covar)[1]

pc_prefix = set_pc_prefix(pheno, covar, ".")
pc_prefix2 = set_pc_prefix(pheno, covar, "_")

# configure logger
Log = Logger(pc_prefix, out_dir2)
Log.just_log("\nvcf2gwas v0.5 \n")
Log.just_log("Initialising..\n")
Log.print_log(f'Start time: {timestamp}\n')
Log.print_log("Parsing arguments..")

#snp_file2 = f'input/{snp_file}'
gene_file = P.set_gene_file()
gene_file_path = None
if gene_file != None:
    gene_file_path = gene_file
    gene_file = os.path.split(gene_file)[1]
gene_thresh = P.set_gene_thresh()

X = P.get_phenotypes()
Y = P.get_covariates()
A = P.set_A()
B = P.set_B()
if A == True and X != None:
    X = None
    Log.print_log("Warning: option 'allphenotypes' will overwrite your phenotype selection")
if B == True and Y != None:
    Y = None
    Log.print_log("Warning: option 'allcovariates' will overwrite your covariate selection")

lm = P.set_lm()
gk = P.set_gk()
eigen = P.set_eigen()
lmm = P.set_lmm()
bslmm = P.set_bslmm()

model = set_model(lm, gk, eigen, lmm, bslmm)
if model != None:
    model2 = model.removeprefix("-")

n_top = P.set_n_top()

# file checking
Log.print_log(f'Genotype file: {snp_file}')
if pheno_files != None:
    Log.print_log(f'Phenotype file(s): {listtostring(pheno_files, ", ")}')
if covar != None:
    Log.print_log(f'Covariate file: {covar}')
if gene_file != None:
    Log.print_log(f'Gene comparison file: {gene_file}')
Log.print_log("Arguments parsed successfully\n")

if model == None:
    Log.print_log("Warning: No model specified for GEMMA analysis")
if A == False and B == False:
    if X == None and Y == None:
        Log.print_log("Warning: No phenotypes (or covariates) specified for GEMMA analysis")

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
check_files3(covar, covar2)

# file preparation
Log.print_log("Preparing files\n")

memory = P.set_memory()
memory_org = memory
memory_total = int(((virtual_memory().total/1e9)//2)*1e3)
if memory_total < memory:
    memory = memory_total
memory2 = memory

threads = P.set_threads()
threads_org = threads
cpu = mp.cpu_count()-1
if cpu == 0:
    cpu = 1
if cpu < threads:
    threads = cpu

n = set_n(lm, gk, lmm, bslmm)
filename = P.set_filename()
if filename != None:
    filename = os.path.split(filename)[1]
min_af = P.set_q()
pca = P.set_pca()
keep = P.set_keep()
multi = P.set_multi()
seed = P.set_seed()

umap_n = P.set_U()
if umap_n != None:
    Starter.check_vals("UMAP", umap_n, 1, 5, Log)
    umap_switch = True
else:
    umap_switch = False

pca_n = P.set_P()
if pca_n != None:
    Starter.check_vals("PCA", pca_n, 2, 10, Log)
    pca_switch = True
else:
    pca_switch = False

pheno_list = []
X_list = []
threads_list = []
l = None
pheno_files2 = []

analysis_num = 0

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
    X = []

#################### Split up files ####################

if pheno_files != None:

    if umap_switch == True and pca_switch == True:
        num = len(pheno_files)*2
    else:
        num = len(pheno_files)

    if num > threads or num > memory/1e3:
        if umap_switch == True and pca_switch == True:
            exit(Log.print_log("Error: Too many input files for available ressources when performing both UMAP and PCA"))
        else:
            exit(Log.print_log("Error: Too many phenotype input files for available ressources"))

    x = 0
    for pheno_file in pheno_files:

        pheno_path = pheno_files_path[x]
        pheno_file2 = os.path.join(pheno_path, pheno_file)
        df = Processing.load_pheno(pheno_file2)
        l = len(df.columns)

        if umap_switch == True:
            if l < umap_n:
                exit(Log.print_log(f'Error: not enough phenotypes in {pheno_file} to perform UMAP (at least {umap_n} required)'))
            else:
                Log.print_log(f"Performing UMAP dimensionality reduction on {pheno_file}..")
                pheno_file3 = pheno_file
                pheno_file3 = pheno_file3.removesuffix(".csv")
                pheno_file3 = f'{pheno_file3}_umap.csv'
                pheno_files2.append(pheno_file3)
                Starter.umap_calc(df, pheno_file3, umap_n, seed, pheno_path)
                Log.print_log(f'Saved as "{pheno_file3}" temporarily in {pheno_path}\nUMAP calculated successful\n')
        
        if pca_switch == True:
            if l < pca_n:
                exit(Log.print_log(f'Error: not enough phenotypes in {pheno_file} to perform PCA (at least {pca_n} required)'))
            else:
                Log.print_log(f"Performing PCA dimensionality reduction on {pheno_file}..")
                pheno_file4 = pheno_file
                pheno_file4 = pheno_file4.removesuffix(".csv")
                pheno_file4 = f'{pheno_file4}_pca.csv'
                pheno_files2.append(pheno_file4)
                Starter.pca_calc(df, pheno_file4, pca_n, pheno_path)
                Log.print_log(f'Saved as "{pheno_file4}" temporarily in {pheno_path}\nPCA calculated successful\n')

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

        if A == False and len(X) > len(df.columns):
            Log.print_log(f"Warning:\nMore phenotypes selected than present in {pheno_file}!\nUsing all available phenotypes..")
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
            Starter.split_phenofile1(X, X_list, df, pheno_file, pheno_list, pheno_path)
            Log.print_log("Phenotype file split up successful")

        else:
            rest = 0
            threads2 = 1
            memory2 = memory//threads
            col_dict = {}
            rest2 = l%threads
            threads3 = l//threads
            Starter.split_phenofile2(threads, threads3, rest2, col_dict, X_list, df, pheno_file, pheno_list, pheno_path)
            Log.print_log("Phenotype file split up successful")
        
        x += 1
        if umap_switch == False and pca_switch == False:
            analysis_num += l

if umap_switch == True and pca_switch == True:
    analysis_num += (umap_n + pca_n)*(len(pheno_files2)//2)

if memory2 < 1000:
    Log.print_log("Warning: Memory might not be sufficient to carry out analysis!")
if memory2 == 0:
    exit(Log.print_log("Error: Memory not sufficient to carry out analysis!"))

# compressing, indexing and filtering VCF file

if snp_file.endswith(".vcf"):
    Log.print_log("Compressing VCF file..")
    timer = time.perf_counter()
    snp_file2 = Converter.compress_snp_file(snp_file2)
    timer_end = time.perf_counter()
    timer_total = round(timer_end - timer, 2)
    snp_file = f'{snp_file}.gz'
    Log.print_log(f'VCF file successfully compressed (Duration: {runtime_format(timer_total)})\n')

try:
    snp_prefix = snp_file.removesuffix(".vcf.gz")
except Exception:
    snp_prefix = snp_file.removesuffix(".gz")

Log.print_log("\nFiltering SNPs..")
snp_file_path = os.path.split(snp_file2)[0]
temp_dir = os.path.join(snp_file_path, "temp")
dir_check = make_dir(temp_dir)
snp_file2_org = snp_file2
snp_file2 = os.path.join(temp_dir, snp_file)

Log.print_log("Indexing VCF file..")
timer = time.perf_counter()
Converter.index_vcf(snp_file2_org)
timer_end = time.perf_counter()
timer_total = round(timer_end - timer, 2)
Log.print_log(f'VCF file successfully indexed (Duration: {runtime_format(timer_total)})')

timer = time.perf_counter()
Converter.filter_snps(min_af, snp_file2_org, snp_file2)
timer_end = time.perf_counter()
timer_total = round(timer_end - timer, 2)
Log.print_log(f'SNPs successfully filtered (Duration: {runtime_format(timer_total)})')
os.remove(f'{snp_file2_org}.csi')

#################### Prepare commands for main.py ####################

args_list = []
args = sys.argv[1:]
# for testing
if args == []:
    args = argvals
args = Starter.delete_string(args, ['-v', '--vcf', '-T', '--threads', '-M', '--memory'])
args.insert(0, snp_file2)
args.insert(0, "--vcf")
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
    Starter.edit_args2(pheno_list, args, args_list, threads_list, pheno, A, X_list, pheno_files_path)
    if umap_switch == True:
        Log.print_log(f'Info:\nAfter reducing dimensions of {pheno} via UMAP, it has been split up in {len(pheno_list)} parts in order to ensure maximum efficiency')
    else:
        Log.print_log(f'Info:\n{pheno} has been split up in {len(pheno_list)} parts in order to ensure maximum efficiency')

else:
    Starter.edit_args3(args, threads, args_list)

Log.print_log("\nFile preparations completed")

#################### Run main.py for analysis ####################

Log.print_log("\nStarting analysis..")
with concurrent.futures.ProcessPoolExecutor(mp_context=mp.get_context('fork'), max_workers=threads_org) as executor:
        executor.map(Starter.run_vcf2gwas, args_list)
Log.print_log("Analysis successfully completed\n")

#################### summary and clean-up ####################

if dir_check == False:
    shutil.rmtree(temp_dir, ignore_errors=True)
else:
    os.remove(snp_file2)
snp_file2 = snp_file2_org

if model == None:
    path = out_dir2
else:
    path = os.path.join(out_dir2, model2)

# summarizer and gene comparison
if model not in ("-gk", "-eigen", None):
    path2 = os.path.join(path, "summary")
    path3 = os.path.join(path2, "top_SNPs")
    make_dir(path2)
    prefix_list = []
    for pheno in pheno_list:
        pc_prefix3 = set_pc_prefix(pheno, covar, "_")
        prefix_list.append(pc_prefix3)
    filenames = Post_analysis.summarizer(path3, path2, pc_prefix3, snp_prefix, n_top, Log, prefix_list)
    if gene_file != None:
        Post_analysis.gene_compare(filenames, gene_file, gene_file_path, gene_thresh, path2, pc_prefix3, snp_prefix, Log)

# move split up files (and reduce "files" folder)
if switch == False and len(pheno_list) > 1:
    for file in pheno_list:
        os.remove(os.path.join(pheno_files_path[0], file))
    if umap_switch == False and pca_switch == False:
        path5 = os.path.join(path, "files")
        pheno_temp = [f'{x.removesuffix(".csv")}_{snp_prefix}' for x in pheno_list]
        for folder in os.listdir(path5):
            if pheno_temp[0] in folder:
                folder2 = folder.replace(".part1", "")
                shutil.rmtree(os.path.join(path5, folder2), ignore_errors=True)
                shutil.move(os.path.join(path5, folder), os.path.join(path5, folder2))
                path6 = os.path.join(path5, folder2)
                for file in os.listdir(path6):
                    for string in [".log.txt", ".cXX.txt", ".log", ".bim", ".bed", ".nosex"]:
                        if file.endswith(string):
                            os.rename(os.path.join(path6, file), os.path.join(path6, file.replace(".part1", "")))
        for folder in os.listdir(path5):
            for pheno in pheno_temp:
                if pheno in folder:
                    for file in os.listdir(os.path.join(path5, folder)):
                        if file.endswith(".fam"):
                            shutil.move(os.path.join(path5, folder, file), os.path.join(path5, folder2, file))
                    shutil.rmtree(os.path.join(path5, folder))

# move umap/pca files
if umap_switch == True or pca_switch ==True:
    switch_names = ["UMAP", "PCA"]
    x = 0
    for switch in [umap_switch, pca_switch]:
        if switch == True:
            path4 = os.path.join(path, switch_names[x])
            make_dir(path4)
            y = 0 
            for pheno_file in pheno_files:
                pheno_path = pheno_files_path[y]
                if umap == True and pca_switch == True:
                    pheno_path = pheno_files_path[y//2]
                for files in os.listdir(pheno_path):
                    if files.startswith(f'{pheno_file.removesuffix(".csv")}_{switch_names[x].lower()}'):
                        shutil.move(os.path.join(pheno_path, files), os.path.join(path4, files))
                y += 1
            Log.print_log(f'Moved {switch_names[x]} files to {path4}')
        x += 1

finish = time.perf_counter()
time_total = round(finish-start, 2)
time_total = runtime_format(time_total)

Log.print_log(f'Clean up successful \n\nvcf2gwas has been successfully completed! \nRuntime: {time_total}\n')

if pheno_files != None:
    pheno_files = listtostring(pheno_files, ', ')
if X != None:
    X1 = []
    if umap_switch == True:
        X1.append(f'{umap_n} UMAP embedding(s)')
    if pca_switch == True:
        X1.append(f'{pca_n} principal components')
    if umap_switch == False and pca_switch == False:
        X1 = X
    X = listtostring(X1, ', ')
if Y != None:
    Y = listtostring(Y, ', ')

Log.summary(snp_file, pheno_files, covar, X, Y, model2, n, filename, min_af, A, B, pca, keep, memory, threads, n_top, gene_file, gene_thresh, multi, umap_n, pca_n, out_dir2, analysis_num)

if model != None:
    shutil.move(os.path.join(out_dir2, f'vcf2gwas{pc_prefix}.log.txt'), os.path.join(out_dir2, model2, f'vcf2gwas_{snp_prefix}{pc_prefix2}.log.txt'))