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

from numpy import log10
from vcf2gwas.parsing import Parser
import warnings

from matplotlib.pyplot import axes, text
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
import subprocess
import sys
import shutil
import logging
import argparse
import multiprocessing as mp
import random
import time

from psutil import virtual_memory
try:
    import numpy as np
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "anaconda", "numpy==1.20*"])
    import numpy as np
try:
    import pandas as pd
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "anaconda", "pandas==1.2*"])
    import pandas as pd
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "conda-forge", "matplotlib==3.4*"])
    import matplotlib.pyplot as plt
try:
    import seaborn as sns
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "anaconda", "seaborn==0.11*"])
    import seaborn as sns
try:
    from sklearn.decomposition import PCA
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "anaconda", "scikit-learn==0.24*"])
    from sklearn.decomposition import PCA
try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "conda-forge", "adjusttext==0.7*"])
    from adjustText import adjust_text
try:
    import umap
except ModuleNotFoundError:
    subprocess.run(["pip", "install", "umap-learn"])
    import umap
try:
    import zarr
except ModuleNotFoundError:
    subprocess.run(["pip", "install", "zarr"])
try:
    import allel
except ModuleNotFoundError:
    subprocess.run(["pip", "install", "scikit-allel"])

from parsing import *

sns.set_style('white')
#sns.set_style('ticks')
#os.chdir(os.path.dirname(os.path.abspath(__file__)))

#set fontsize
argvals = None
P = Parser(argvals)
fontsize = P.set_fontsize()
fontsize2 = fontsize - fontsize*0.1
fontsize3 = fontsize - fontsize*0.5
fontsize4 = fontsize - fontsize*0.7
fontsize5 = fontsize - fontsize*0.92
fontsize6 = fontsize - fontsize*0.94

plt.rc("lines", linewidth=fontsize5, markersize=fontsize4)
plt.rc('font', size=fontsize3, weight="bold")
plt.rc('axes', linewidth=fontsize5)
plt.rc('xtick.major', width=fontsize5, size=fontsize4)
plt.rc('ytick.major', width=fontsize5, size=fontsize4)
sns.set_color_codes()

############################## Functions ##############################

def flatten_list(l):

    l = [item for sublist in l for item in sublist]
    l = list(dict.fromkeys(l))
    return l

def set_pc_prefix(pheno_file, covar_file, string):
    """Description:
    set variable to prefix from pheno/covar input file"""

    if covar_file != None:
        pc_prefix = f'{string}{covar_file.removesuffix(".csv")}'
    if pheno_file != None:
        pc_prefix = f'{string}{pheno_file.removesuffix(".csv")}'
    if covar_file == None and pheno_file == None:
        pc_prefix = string
    return pc_prefix

def write_timestamp(time):

    file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_timestamp.txt"), 'a')
    file.write(f'{time}\n')
    file.close()

def read_timestamp():

    file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_timestamp.txt"), "r")
    lines = file.readlines()
    last = lines[-1].rstrip()
    return last      

def check_files(snp_file2, gene_file, gene_file_path):
    """Description:
    checks if input files exists (SNP file and Gene file)"""

    if os.path.isfile(snp_file2) == False:
        sys.exit(print("Error: VCF file non existent"))
    if gene_file != None:
        if os.path.isfile(gene_file_path) == False: sys.exit(print("Error: Gene comparison file non existent, please check if file path was specified correctly"))

def check_files2(pheno_file, pheno_file2):
    """Description:
    checks if input files exists (Pheno file)"""

    if pheno_file != None:
        if os.path.isfile(pheno_file2) == False: sys.exit(print("Error: Phenotype file non existent, please check if file path was specified correctly"))

def check_files3(covar_file, covar_file2):
    """Description:
    checks if input files exists (Covar file)"""

    if covar_file != None:
        if os.path.isfile(covar_file2) == False: sys.exit(print("Error: Covariate file non existent, please check if file path was specified correctly"))

def make_dir(name):
    """Description: 
    creates the directory "name" if it doesn't already exist"""

    try:
        dir_check = False
        os.makedirs(name)
    except OSError as e:
        dir_check = True
    return dir_check

def listtostring(l, string=' '):
    """Description:
    converts list to string"""

    try:
        return string.join(map(str, l))
    except Exception as e:
        pass

def concat_lists(X, Y):
    """Description:
    concatenates length of two lists and returns as list of consecutive integers"""

    c = len(X) + len(Y)
    return list(range(1, c+1))

def set_cols(cols, X, pheno_subset):
    """Description:
    returns list of column names"""
    
    if pheno_subset.empty == True:
        return []
    else:
        for i in X:
            name = pheno_subset.columns[i-1]
            cols.append(name)
        return cols

def runtime_format(time_total):
    """Description:
    returns runtime formatted to hour, minute and second format"""

    h = m = s = ""
    if time_total >= 3600:
        n = time_total//3600
        time_total = time_total - 3600*n
        h = f'{int(n)} hours, '
        if n == 1:
            h = f'{int(n)} hour, '
    if time_total >= 60:
        n = time_total//60
        time_total = time_total - 60*n
        m = f'{int(n)} minutes, '
        if n == 1:
            m = f'{int(n)} minute, '
    time_total = round(time_total, 1)
    s = f'{time_total} seconds'
    if time_total == 1.0:
        s = f'{time_total} second'
    return h+m+s

def pheno_switcher(genos, phenos):

    phenos_org = phenos
    try:
        phenos = [int(i) for i in phenos]
    except:
        pass

    for g in genos:

        df = pd.read_csv(g, index_col=0)
        x = df.columns
        try:
            if any(i > len(x) for i in phenos):
                phenos = phenos_org
        except:
            pass
        x_int = [i for i in range(1, len(x)+1)]
        pdict = dict(zip(x, x_int))

        if all(i in x for i in phenos):
            switch = 1
        elif any(i in x for i in phenos):
            switch = -1
        else:
            if all(i in pdict.values() for i in phenos):
                switch = 0
            else:
                switch = -1

        if switch == -1:
            print("Error: Please make sure to only select existing phenotypes / covariates")
            sys.exit(1)
        elif switch == 1:
            phenos = [pdict[i] for i in phenos]
        elif switch == 0:
            pass

    return phenos

def pheno_switcher2(genos, phenos, phenos2):

    phenos_org = phenos
    try:
        phenos = [int(i) for i in phenos]
    except:
        pass

    if phenos == None:
        phenos = []

    for g in genos:

        df = pd.read_csv(g, index_col=0)
        x = df.columns
        x_list = list(x.values)
        x_int = [i for i in range(1, len(x)+1)]
        try:
            if any(i > len(x) for i in phenos):
                phenos = phenos_org
        except:
            pass
        pdict = dict(zip(x, x_int))
        pdict2 = dict(zip(x_int, x))

        if all(i in x for i in phenos):
            switch = 1
        elif any(i in x for i in phenos):
            switch = -1
        elif phenos == None:
            switch = 1
        else:
            if all(i in pdict.values() for i in phenos):
                switch = 0
            else:
                switch = -1    

        if switch == -1:
            print("Error: Please make sure to only select existing phenotypes / covariates")
            sys.exit(1)
        elif switch == 1:
            phenos2 = [pdict2[i] for i in phenos2]
            x_list = phenos2
        elif switch == 0:
            x_list = [pdict2[i] for i in phenos2]

    return phenos2, x_list

def set_model(lm, gk, eigen, lmm, bslmm):

    l = ["-lm", "-gk", "-eigen", "-lmm", "-bslmm"]
    model = None
    n = 0
    for i in (lm, gk, eigen, lmm, bslmm):
        if i != None and i != False:
            model = l[n]
        n += 1
    return model

def set_n(lm, gk, lmm, bslmm):

    n = None
    for i in (lm, gk, lmm, bslmm):
        if i != None:
            n = str(i)
    return n

def move_log(model, model2, pc_prefix, snp_prefix, out_dir, timestamp):

    if model != None:
        path = os.path.join(out_dir, model2, "logs", "temp")
        make_dir(path)
        if model in ["-gk", "-eigen"]:
            shutil.move(os.path.join(out_dir, f'vcf2gwas{pc_prefix}.log.txt'), os.path.join(path, f'vcf2gwas_analysis_{snp_prefix}{pc_prefix}_{timestamp}.log.txt'))
        else:
            shutil.move(os.path.join(out_dir, f'vcf2gwas{pc_prefix}.log.txt'), os.path.join(path, f'vcf2gwas_analysis_{snp_prefix}_{pc_prefix}_{timestamp}.log.txt'))


############################## Classes ##############################

class Logger:
    """Description:
    contains functions regarding the log and output messages"""

    def __init__(self, pc_prefix, out_dir):
        """Description:
        initializes logger"""

        LOG_FORMAT = "%(message)s"
        logging.basicConfig(filename=os.path.join(out_dir, f'vcf2gwas{pc_prefix}.log.txt'), level=logging.INFO, format=LOG_FORMAT, filemode='w')
        self.logger = logging.getLogger()

    def print_log(self, text):
        """Description:
        prints text to both log file and command-line"""

        print(text)
        self.logger.info(text)

    def just_log(self, text):
        """Description:
        prints text to log file"""        

        self.logger.info(text)
    
    def summary(
        self, snp_file, pheno_file, covar_file, X, Y, model, N, filename, min_af, A, B, 
        pca, keep, memory, threads, n_top, gene_file, species, gene_thresh, multi, umap_n, pca_n, 
        out_dir, analysis_num, sigval, nolabel, chr, chr2, chr_num, X_names, snp_total, snp_sig, sig_level, geno_pca_switch, burn, sampling, snpmax, noqc,
        input_str
    ):
        """Description:
        prints summary of input variables and methods"""

        a = b = c = d = e = f = g = h = i = j = k = l = m = n = o = p = q = r = s = t = u = v = w = x = y = z = aa = ab = ac = ad = ae = af = ag = ah = aj = ak = al = ""

        model_dict = {
            "lm" : "linear model",
            "gk" : "estimate relatedness matrix",
            "eigen" : "eigen-decomposition of relatedness matrix",
            "lmm" : "linear mixed model",
            "bslmm" : "bayesian sparse linear mixed model"
        }
        lmm_method_dict = {
            "1" : "Wald test",
            "2" : "likelihood ratio test",
            "3" : "score test",
            "4" : "Wald test, likelihood ratio test, score test"
        }
        bslmm_method_dict = {
            "1" : "standard linear BSLMM",
            "2" : "ridge regression/GBLUP",
            "3" : "probit BSLMM"
        }
        gk_method_dict = {
            "1" : "centered relatedness matrix",
            "2" : "standardized relatedness matrix"
        }


        al = f'\n{input_str}'
        v = f'\n{out_dir}'
        w = f'\n{analysis_num}'
        if X_names != "":
            ab = f'({X_names})'
        z = f'\n{chr_num}'
        ac = chr2
        ad = snp_total
        ae = snp_sig
        af = sig_level
        a = f'\n- VCF file: "{snp_file}"'
        if pheno_file != None:
            b = f'\n- Phenotype file(s): "{pheno_file}"'
            if X == "":
                c = '\n  Phenotypes chosen: None'
            else:
                c = f'\n  Phenotypes chosen: {X}'
        if covar_file != None:
            d = f'\n- Covariate file: "{covar_file}"' 
            if Y == "":
                e = '\n  Covariates chosen: None'
            else:
                e = f'\n  Covariates chosen: {Y}'
            if geno_pca_switch == True:
                d = f'\n- Covariates: principal components ({Y})'
                e = "" 
        if gene_file != None:
            if species == None:
                q = f'\n- Gene comparison file: "{gene_file}"'
            else:
                q = f'\n- Gene comparison with: {species}'
            if gene_thresh != 100000:
                r = f'\n- Gene comparison threshold: {gene_thresh}'
        if model != None:
            f = f'\n- Model: "{model}" ({model_dict[model]})'
        if N != None:
            if model in ["lm", "lmm"]:
                g = f'\n- Test method: "{N}" ({lmm_method_dict[N]}) '
            if model == "gk":
                g = f'\n- Test method: "{N}" ({gk_method_dict[N]}) '
            if model == "bslmm":
                g = f'\n- Test method: "{N}" ({bslmm_method_dict[N]}) '
        if filename != None:
            h = f'\n- Relatedness matrix filename: "{filename}"'
        if min_af != 0.01:
            i = f'\n  --minaf {min_af}'
        if chr != None:
            aa = f'\n  --chromosome {chr2}'
        if n_top != 15:
            j = f'\n  --topsnp {n_top}'
        if memory != int(((virtual_memory().total/1e9)//2)*1e3):
            k = f'\n  --memory {memory}'
        if threads != mp.cpu_count()-1:
            l = f'\n  --threads {threads}'
        if sigval != None:
            x = f'\n  --sigval {sigval}'
        if burn != 100000:
            ag = f'\n  --burn {burn}'
        if sampling != 1000000:
            ah = f'\n  --sampling {sampling}'
        if snpmax != 300:
            aj = f'\n  --snpmax {snpmax}'
        if A == True:
            m = '\n  --allphenotypes'
            c = f'\n  Phenotypes chosen: all phenotypes'
        if B == True:
            if geno_pca_switch == False:
                n = '\n  --allcovariates'
                e = f'\n  Covariates chosen: all covariates'
        if keep == True:
            o = '\n  --retain'
        if pca != None:
            p = f'\n  --kcpca (threshold: {pca})'
        if multi == True:
            s = '\n  --multi'
            f = f'\n- Model: "{model}" (multivariate linear mixed model)'
        if nolabel == True:
            y = '\n  --nolabel'
        if noqc == True:
            ak = '\n  --noqc'
        if umap_n != None:
            t = '\n  --UMAP'
            if umap_n != 2:
                t = f'\n  --UMAP {umap_n}'
        if pca_n != None:
            u = '\n  --PCA'
            if pca_n != 2:
                u = f'\n  --PCA {pca_n}'

        self.logger.info(
            f'\nSummary:\n\nOutput directory:{v}\n\nPhenotypes analyzed in total:{w} {ab}\n\nChromosomes analyzed in total:{z} ({ac})\n\nVariants analyzed: \nTotal: {ad} \nSignificant: {ae} \nLevel of significance: {af} \n\n\nInput:\n\nCommand:{al}\n\nFiles:{a}{b}{c}{d}{e}{q}{r}{h}\n\nGEMMA parameters:{f}{g}\n\nOptions:{t}{u}{s}{i}{aa}{j}{k}{l}{x}{ag}{ah}{aj}{m}{n}{y}{ak}{o}{p}'
        )


class Starter:
    """Description:
    contains functions regarding the starter script taking care of multiprocessing"""

    def run_vcf2gwas(args):
        process = subprocess.run(args)
        code = process.returncode
        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"), 'a')
        file.write(str(code))
        file.close()

    def get_gene_file(gene_file):

        gene_dict = {
            "HS" : "Homo sapiens",
            "MM" : "Mus musculus",
            "DR" : "Danio rerio",
            "RN" : "Rattus norvegicus",
            "GG" : "Gallus gallus",
            "DM" : "drosophila melanogaster",
            "CE" : "Caenorhabditis elegans",
            "SC" : "Saccharomyces cerevisiae",
            "AG" : "Anopheles gambiae",
            "AT" : "Arabidopsis thaliana",
            "SL" : "Solanum lycopersicum",
            "ZM" : "Zea mays",
            "OS" : "Oryza sativa",
            "VV" : "Vitis vinifera"
        }

        gene_file_dict = {
            "HS" : "GFF_HS.csv",
            "MM" : "GFF_MM.csv",
            "DR" : "GFF_DR.csv",
            "RN" : "GFF_RN.csv",
            "GG" : "GFF_GG.csv",
            "DM" : "GFF_DM.csv",
            "CE" : "GFF_CE.csv",
            "SC" : "GFF_SC.csv",
            "AG" : "GFF_AG.csv",
            "AT" : "GFF_AT.csv",
            "SL" : "GFF_SL.csv",
            "ZM" : "GFF_ZM.csv",
            "OS" : "GFF_OS.csv",
            "VV" : "GFF_VV.csv"
        }

        if gene_file in gene_dict:
            species = gene_dict[gene_file]
            gene_file = os.path.join(os.path.dirname(__file__), "GFF_files", gene_file_dict[gene_file])
        else:
            species = None

        return gene_file, species

    def check_vals(str, var, min, max, Log):
        if var == 0:
            sys.exit(Log.print_log(f"Error: input value for dimension reduction ({str}) must be greater than zero!"))
        elif var < min:
            Log.print_log(f'Warning: input value for dimension reduction ({str}) is not recommended to be less than {min}!\n')
        elif var > max:
            Log.print_log(f'Warning: input value for dimension reduction ({str}) is not recommended to be greater than {max}!\n')

    def adjust_threads(list, threads, rest, threads_list):
        for i in range(len(list)):
            if i+1 <= rest:
                threads2 = threads + 1
            else:
                threads2 = threads
            threads_list.append(str(threads2))

    def label_point(x, y, val, texts):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            texts.append(plt.text(point['x'], point['y']+np.random.random()/100, str(point['val']), fontsize="medium", fontweight="bold"))

    def umap_calc(df, pheno_file, n, seed, pheno_path, Log):
        """Description:
        performs dimension reduction via UMAP"""

        df.fillna("NA", inplace=True)
        df = df.replace("NA",np.nan)
        df = df.dropna(how="any")
        if df.empty:
            sys.exit(Log.print_log("Error: Too many missing values to perform UMAP!"))
        x = df.values
        y = df.index.tolist()
        if seed == True:
            rand = random.randint(1, 1000)
        else:
            rand = 322
        emb_list = []
        for i in range(n):
            emb_list.append(f'Emb{i+1}')
        
        #run umap 
        reducer = umap.UMAP(random_state=rand, min_dist=0.1, n_components=n, n_neighbors=5, metric="manhattan")
        embedding = reducer.fit_transform(x)
        embeddingDf = pd.DataFrame(data = embedding, columns=emb_list )
        embeddingDf['host']=y

        #plot embeddings
        if n == 2:
            texts = []
            plt.subplots(figsize=(12,12))
            sns.scatterplot("Emb1","Emb2", data=embeddingDf)
            sns.despine(offset=10)
            Starter.label_point(embeddingDf.Emb1, embeddingDf.Emb2, embeddingDf.host, texts) #custom function used
            plt.title("UMAP embeddings", fontsize=fontsize)
            plt.xticks(fontsize=fontsize2, fontweight="bold")
            plt.yticks(fontsize=fontsize2, fontweight="bold")
            plt.xlabel("Embedding 1", fontsize=fontsize, fontweight="bold")
            plt.ylabel("Embedding 2", fontsize=fontsize, fontweight="bold")
            adjust_text(texts, embeddingDf['Emb1'].values, embeddingDf['Emb2'].values, expand_text=(1.02, 1.02), expand_align=(1.02, 1.02), force_text=(0,0.7), lim=250, arrowprops=dict(arrowstyle="-", color='k', lw=0.5, alpha=0.6))
            plt.savefig(os.path.join(pheno_path, f'{pheno_file.removesuffix(".csv")}.png'))
            plt.close()

        embeddingDf.set_index("host", inplace=True, drop=True)
        embeddingDf.index.name = None
        embeddingDf.to_csv(os.path.join(pheno_path, pheno_file))

    def pca_calc(df, pheno_file, n, pheno_path, Log):
        """Description:
        performs principal component analysis"""

        df.fillna("NA", inplace=True)
        df = df.replace("NA",np.nan)
        df = df.dropna(how="any")
        if df.empty:
            sys.exit(Log.print_log("Error: Too many missing values to perform PCA!"))
        x = df.values
        y = df.index.tolist()

        pc_list = []
        for i in range(n):
            pc_list.append(f'PC{i+1}')

        #run pca 
        pca = PCA(n_components=n)
        embedding = pca.fit_transform(x)
        embeddingDf = pd.DataFrame(data = embedding, columns=pc_list)
        embeddingDf['host']=y
        
        #plot scree plot
        df_scree = pd.DataFrame({'var':pca.explained_variance_ratio_, 'PC':pc_list})
        plt.subplots(figsize=(16,12))
        sns.barplot(x='PC', y='var', data=df_scree, color='b')
        plt.title("Scree plot", fontsize=fontsize)
        plt.xticks(fontsize=fontsize2, fontweight="bold")
        plt.yticks(fontsize=fontsize2, fontweight="bold")
        plt.xlabel("Principal Components", fontsize=fontsize, fontweight="bold")
        plt.ylabel("Variance Explained", fontsize=fontsize, fontweight="bold")
        plt.savefig(os.path.join(pheno_path, f'{pheno_file.removesuffix(".csv")}_scree.png'))
        plt.close()

        #plot PCs
        if n == 2:
            texts = []
            plt.subplots(figsize=(12,12))
            sns.scatterplot("PC1","PC2", data=embeddingDf)
            sns.despine(offset=10)
            Starter.label_point(embeddingDf.PC1, embeddingDf.PC2, embeddingDf.host, texts) #custom function used
            plt.title("Principal components", fontsize=fontsize)
            plt.xticks(fontsize=fontsize2, fontweight="bold")
            plt.yticks(fontsize=fontsize2, fontweight="bold")
            plt.xlabel("Principal Component 1", fontsize=fontsize, fontweight="bold")
            plt.ylabel("Principal Component 2", fontsize=fontsize, fontweight="bold")
            adjust_text(texts, embeddingDf['PC1'].values, embeddingDf['PC2'].values, expand_text=(1.02, 1.02), expand_align=(1.02, 1.02), force_text=(0,0.7), lim=250, arrowprops=dict(arrowstyle="-", color='k', lw=0.5, alpha=0.6))
            plt.savefig(os.path.join(pheno_path, f'{pheno_file.removesuffix(".csv")}.png'))
            plt.close()

        embeddingDf.set_index("host", inplace=True, drop=True)
        embeddingDf.index.name = None
        embeddingDf.to_csv(os.path.join(pheno_path, pheno_file))

    def split_phenofile1(X, X_list, df, pheno_file, pheno_list, pheno_path):

        ind = 1
        for i in X:
            df_new = df.iloc[:,i-1]
            df_new.fillna("NA", inplace=True)
            df_name = pheno_file.removesuffix(".csv")
            df_name = f'{df_name}.part{ind}.csv'
            pheno_list.append(df_name)
            df_new.to_csv(os.path.join(pheno_path, df_name))
            ind += 1

    def split_phenofile2(threads, threads3, rest2, col_dict, X_list, df, pheno_file, pheno_list, pheno_path):

        for i in range(threads):
            x = 0
            cols = []
            for c in range(threads3):
                cols.append(i+x)
                x = x+threads
            if i+1 <= rest2:
                cols.append(i+x)
            col_dict[i] = cols

            df_new = df.iloc[:,col_dict[i]]
            df_new.fillna("NA", inplace=True)
            df_name = pheno_file.removesuffix(".csv")
            df_name = f'{df_name}.part{i+1}.csv'
            pheno_list.append(df_name)
            df_new.to_csv(os.path.join(pheno_path, df_name))

    def delete_string(args, strings):

        for string in strings:
            temp_list = [i for i, e in enumerate(args) if e == string]
            temp_list2 = [x+1 for x in temp_list]
            temp_list = temp_list+temp_list2
            args = [j for i, j in enumerate(args) if i not in temp_list]
        return args

    def edit_args1(pheno_list, args, args_list, threads_list, umap_switch, pca_switch, A, pheno_path):

        for i in range(len(pheno_list)):
            args2 = list(args)
            args_list.append(args2)
            args_list[i].insert(4, '--threads')
            args_list[i].insert(5, threads_list[i])
            args_list[i].insert(6, "--pfile")
            if umap_switch == True and pca_switch == True:
                args_list[i].insert(7, os.path.join(pheno_path[i//2], pheno_list[i]))
            else:
                args_list[i].insert(7, os.path.join(pheno_path[i], pheno_list[i]))
            if umap_switch == True or pca_switch == True:
                if A == False:
                    args_list[i].insert(8, "--allphenotypes")

    def edit_args2(pheno_list, args, args_list, threads_list, pheno_file, A, X_list, pheno_path):
        
        for i in range(len(pheno_list)):
            args2 = list(args)
            args_list.append(args2)
            args_list[i].insert(4, '--threads')
            args_list[i].insert(5, threads_list[i])
            args_list[i] = [os.path.join(pheno_path[0], pheno_list[i]) if pheno_file in x else x for x in args_list[i]]
            #args_list[i] = [os.path.join(pheno_path[0], pheno_list[i]) if x==os.path.join(pheno_path[0], pheno_file) else x for x in args_list[i]]
            #n = 7
            if A == False:
                args_list[i].insert(6, "--allphenotypes")

    def edit_args3(args, threads, args_list):

        args.insert(4, '--threads')
        args.insert(5, str(threads))
        args_list.append(args)

    def check_return_codes(Log):

        os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_timestamp.txt"))
        code_file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"), 'r')
        code_file_str = code_file.read()
        if "1" in code_file_str:
            code_file.close()
            shutil.rmtree("_vcf2gwas_temp", ignore_errors=True)
            #os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"))
            sys.exit()
        if "0" not in code_file_str:
            code_file.close()
            shutil.rmtree("_vcf2gwas_temp", ignore_errors=True)
            #os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"))
            sys.exit(Log.print_log("Error: During the analysis, vcf2gwas encountered an unexpected error"))
        else:
            code_file.close()
            os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"))
            Log.print_log("Analysis successfully completed\n")
            
    def get_snpcounts():

        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_snpcount_total.txt"), "r")
        file2 = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_snpcount_sig.txt"), "r")
        try:
            lines = file.readlines()
            lines2 = file2.readlines()
            file.close()
            file2.close()
            lines = [int(i.rstrip()) for i in lines]
            lines2 = [int(i.rstrip()) for i in lines2]
            x1 = min(lines)
            x2 = max(lines)
            n1 = min(lines2)
            n2 = max(lines2)
        except:
            x1 = x2 = 0
            n1 = n2 = 0
        if x1 == x2:
            x = x1
        else:
            x = listtostring([x1, x2], " - ")
        if n1 == n2:
            n = n1
        else:
            n = listtostring([n1, n2], " - ")
        os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_snpcount_total.txt"))
        os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_snpcount_sig.txt"))
        return x, n 

    def get_sig_level():

        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_sig_level.txt"), "r")
        try:
            lines = file.readlines()
            file.close()
            lines = [i.rstrip() for i in lines]
            x1 = min(lines)
            x2 = max(lines)
        except:
            x1 = x2 = 0
        if x1 == x2:
            x = x1
            if x == 0:
                x = "-"
        else:
            x = listtostring([x1, x2], " - ")
        os.remove(os.path.join("_vcf2gwas_temp", "vcf2gwas_sig_level.txt"))
        return x


class QC:

    def pheno_QC(df, X, folder):

        for x in X:

            x_name = df.columns[x-1]

            fig = plt.figure(figsize=(16,12))
            sns.histplot(df[x_name])
            sns.despine(offset=10)
            plt.title(f'Distribution of {x_name}', fontsize=fontsize)
            plt.xticks(fontsize=fontsize2, fontweight="bold")
            plt.xlabel(x_name, fontsize=fontsize, fontweight="bold")
            plt.yticks(fontsize=fontsize2, fontweight="bold")
            plt.ylabel("Count", fontsize=fontsize, fontweight="bold")
            plt.savefig(os.path.join(folder, f'Distribution_{x_name}.png'))
            plt.close()

    def plot_windowed_variant_density(pos, window_size, outname, qc_dir, title=None):

        # setup windows 
        bins = np.arange(0, pos.max(), window_size)
        
        # use window midpoints as x coordinate
        x = (bins[1:] + bins[:-1])/2
        
        # compute variant density in each window
        h, _ = np.histogram(pos, bins=bins)
        y = h / window_size
        
        # plot
        plt.subplots(figsize=(16, 10))
        sns.lineplot(x=x, y=y)
        sns.despine(offset=10)
        #sns.set_context("talk")
        if title:
            plt.title(title, fontsize=fontsize)
        plt.xticks(fontsize=fontsize2, fontweight="bold")
        plt.yticks(fontsize=fontsize2, fontweight="bold")
        
        plt.xlabel('Chromosome position (bp)', fontsize=fontsize, fontweight="bold")
        plt.ylabel('Variant density (bp$^{-1}$)', fontsize=fontsize, fontweight="bold")
        plt.savefig(os.path.join(qc_dir, f'{outname}.png'))
        plt.close()

    def plot_variant_hist(variants, f, outname, qc_dir, bins=30):
        
        if f == 'DP':
            bins = 50
        x = variants[f][:]
        plt.subplots(figsize=(16, 12))
        sns.histplot(data=x, bins=bins)
        sns.despine(offset=10)
        plt.xticks(fontsize=fontsize2, fontweight="bold")
        plt.yticks(fontsize=fontsize2, fontweight="bold")
        plt.xlabel(f, fontsize=fontsize, fontweight="bold")
        plt.ylabel('No. variants', fontsize=fontsize, fontweight="bold")
        plt.title(f'Variant {f} distribution', fontsize=fontsize)
        plt.savefig(os.path.join(qc_dir, f'{outname}.png'))
        plt.close()

    def geno_QC(vcf, dir_name, qc_dir, chr, Log):

        for i in chr:
            Log.print_log(f'QC for Chromosome: {i}')
            chr_dir = os.path.join(qc_dir, i)
            os.makedirs(chr_dir)
            zname = os.path.join(dir_name, f'vcf_temp_{i}.zarr')
            allel.vcf_to_zarr(vcf, zname, fields='*', region=i, overwrite=True)
            callset = zarr.open_group(zname, mode='r')

            names = ['POS', 'REF', 'ALT']
            for s in ['DP', 'MQ', 'QD']:
                try:
                    callset[f'variants/{s}']
                    names.append(s)
                except Exception:
                    pass

            variants = allel.VariantChunkedTable(callset['variants'], names=names)

            pos = variants['POS'][:]
            QC.plot_windowed_variant_density(pos, 100000, f'VD_{i}', chr_dir, title='Raw variant density')
    
            for val in ['DP', 'MQ', 'QD']:
                if val in names:
                    QC.plot_variant_hist(variants, val, f'{val}_{i}', chr_dir)
            
            shutil.rmtree(zname, ignore_errors=True)


class Processing:
    """Description:
    contains functions regarding processing and editing genotypes and phenotypes""" 

    def load_pheno(pheno_file):
        """Description:
        reads phenotype file and removes duplicates"""

        df = pd.read_csv(pheno_file, index_col=0)
        df = df[~df.index.duplicated(keep="first")]
        return df

    def process_snp_file(snp_file):
        """Description:
        reads accession names from VCF file and returns them as list"""

        out = subprocess.run(['bcftools', 'query', '-l', snp_file], stdout=subprocess.PIPE, text=True)
        ls = pd.Series(out.stdout)
        ls = ls.str.split(pat="\n", expand=True).T
        ls = ls[:-1].astype(str)
        ls = ls.set_index([0])
        ls = ls.rename_axis('acc')
        ls.index = ls.index.str.strip()
        ls = ls.reset_index()
        return ls['acc'].tolist()

    def rm_geno(diff, subset, snp_file):
        """Description:
        removes excess individuals from VCF file"""

        df2 = pd.Series(diff)
        df2.to_csv(f'{subset}.txt', header=None, index=None)

        if df2.empty == True:
            shutil.copy(snp_file, f'{subset}.vcf.gz')
        else:
            remove = subprocess.run(['bcftools', 'view', '-S', f'^{subset}.txt', snp_file, '-Oz', '-o', f'{subset}.vcf.gz'], stdout=subprocess.PIPE, text=True)
            remove.stdout

    def rm_pheno(df, diff, File):
        """Description:
        removes excess individuals from csv file"""

        pheno2 = df[~df.index.isin(diff)]
        #pheno2.to_csv(f'sub_{File}')
        return pheno2

    def pheno_index(file):
        """Description:
        returns list of dataframe index"""

        file.index = file.index.astype(str).str.strip()
        return file.index.tolist()

    def make_diff(l1, l2):
        """Description:
        returns difference of two sets as list"""

        return list(set(l1)-set(l2))

    def make_uniform(l1, l2, diff1, diff2, df, subset, snp_file, File, string, Log):
        """Description:
        removes excess individuals from both VCF and csv file, returns adjusted phenotype data as dataframe"""
        
        if len(l2) == len(l1):
            if set(l2) != set(l1):
                Log.print_log(f'Not all individuals in {string} and genotype file match')
                Processing.rm_geno(diff1, subset, snp_file)
                pheno_subset = Processing.rm_pheno(df, diff2, File)
            else:
                Log.print_log(f'All {string} and genotype individuals match')
                Processing.rm_geno(diff1, subset, snp_file)
                pheno_subset = Processing.rm_pheno(df, diff2, File)
        else:
            Log.print_log(f'Not all individuals in {string} and genotype file match')
            Processing.rm_geno(diff1, subset, snp_file)
            pheno_subset = Processing.rm_pheno(df, diff2, File)
        return pheno_subset

    def prepare_fam(subset2):
        """Description:
        prepares .fam file for editing, removes 5th column"""

        fam = pd.read_csv(f'{subset2}.fam', index_col=0, header=None, sep="\s+")
        fam.index = fam.index.astype(str).str.strip()
        return fam.drop(fam.columns[4], axis=1)

    def edit_fam(fam, pheno_subset, subset2, X, char, string, Log, model, model2, pc_prefix):
        """Description:
        sorts phenotype data according to .fam file then adds selected data to .fam file"""

        pheno_subset_new = pheno_subset.reindex(fam.index)
        for i in X:
            try:
                fam[f'{i}{char}'] = pheno_subset_new.iloc[:,(i-1)]
            except IndexError:
                for file in os.listdir():
                    if file.startswith(subset2):
                        os.remove(file)
                move_log(model, model2, pc_prefix)
                Log.print_log(f'Warning: selected {string}(s) not in {string} file')
                sys.exit(1)
        fam.fillna("NA", inplace=True)
        fam2 = fam
        fam2.dropna(axis=1, how="all", inplace=True)
        if len(fam2.columns) <= 4:
            Log.print_log("Error: Please make sure that the individual names don't contain any special characters/delimiters/spaces")
            sys.exit(1)
        fam.to_csv(f'{subset2}.fam', sep=' ', header=False)
        Log.print_log(f'{string}(s) added to .fam file')
        return pheno_subset_new.columns.tolist()

    def make_covarfile(fam, pheno_subset, subset2, Y, Log):
        """Description:
        creates covariate file for linear mixed model by adding intercept column and removing old index"""

        pheno_subset_new = pheno_subset.reindex(fam.index)
        new = pd.DataFrame()
        for i in Y:
            try:
                new[str(i)] = pheno_subset_new.iloc[:,(i-1)]
            except IndexError:
                for file in os.listdir():
                    if file.startswith(subset2):
                        os.remove(file)
                Log.print_log("Error: selected covariate(s) not in covariate file")
                sys.exit(1)
        new["newindex"] = 1
        new.set_index("newindex", inplace=True)
        new.fillna("NA", inplace=True)
        name = f'{subset2}.covariates.txt'
        new.to_csv(name, sep=" ", header=False)
        return name

    def pca_analysis(subset2, pca, memory, threads, chrom):
        """Description:
        performs LD pruning and principal component analysis via plink"""

        fam = pd.read_csv(f'{subset2}.fam', index_col=0, header=None, sep="\s+")
        n = fam.shape[0]
        del fam
        
        for i in range(5):
            print(f'\nPruning iteration {str(i)}..')
            subprocess.run(['plink', '--bfile', subset2, '--mind', '1', '--indep-pairwise', '100', '10', str(pca), '--allow-no-sex', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)])
            subprocess.run(['plink', '--bfile', subset2, '--mind', '1', '--extract', 'plink.prune.in', '--make-bed', '--out', subset2,'--allow-no-sex', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)])

        print("\nExtracting principal components..")
        do_pca = subprocess.run(['plink', '--pca', str(n), '--bfile', subset2, '--mind', '1', '--out', subset2, '--allow-no-sex', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)], stdout=subprocess.PIPE, text=True)
        do_pca.stdout

        df = pd.read_csv(f'{subset2}.eigenvec', header=None, sep="\s+")
        df = df.drop([0,1], axis=1)

        df.to_csv(f'{subset2}.eigenvec', index=False, header=False,  sep=" " )

        for file in os.listdir():
            if file.endswith("~"):
                os.remove(file)
            if file.startswith("plink"):
                os.remove(file)

    def pca_analysis2(snp_file, n, memory, threads, chrom, list1, dir_temp):
        """Description:
        performs LD pruning and principal component analysis via plink"""
        
        name = os.path.join(dir_temp, "vcf2gwas_geno_pca")
        string = "_"
        list2 = [l for l in list1 if string in l]
        if list2 != None:
            if chrom <= 24:
                do_pca = subprocess.run(['plink', '--pca', str(n), '--vcf', snp_file, '--mind', '1', '--out', name, '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--double-id', '--memory', str(memory),'--threads', str(threads)], stdout=subprocess.PIPE, text=True)
                do_pca.stdout
            else:
                do_pca = subprocess.run(['plink', '--pca', str(n), '--vcf', snp_file, '--mind', '1', '--out', name, '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--double-id', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)], stdout=subprocess.PIPE, text=True)
                do_pca.stdout
        else:
            if chrom <= 24:
                do_pca = subprocess.run(['plink', '--pca', str(n), '--vcf', snp_file, '--mind', '1', '--out', name, '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--memory', str(memory),'--threads', str(threads)], stdout=subprocess.PIPE, text=True)
                do_pca.stdout
            else:
                do_pca = subprocess.run(['plink', '--pca', str(n), '--vcf', snp_file, '--mind', '1', '--out', name, '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)], stdout=subprocess.PIPE, text=True)
                do_pca.stdout        

        df = pd.read_csv(f'{name}.eigenvec', header=None, sep="\s+")
        while len(df.columns) > n+1:
            df = df.drop(df.columns[0], axis=1)
        cols = [""]
        cols2 = []
        for i in range(n):
            cols.append(f'PC{i+1}')
            cols2.append(f'PC{i+1}')
        df.columns = cols
        for file in os.listdir(dir_temp):
            if file.startswith("vcf2gwas_geno_pca"):
                os.remove(os.path.join(dir_temp, file))
        df.to_csv(f'{name}.csv', index=False)
        return cols2


class Converter:
    """Description:
    contains functions regarding formatting and converting input files"""

    def compress_snp_file(snp_file):
        """Description:
        compresses VCF file"""

        subprocess.run(['bcftools', 'view', snp_file, '-Oz', '-o', f'{snp_file}.gz'])
        return f'{snp_file}.gz'

    def index_vcf(snp_file):
        """Description:
        indexes VCF file"""

        process = subprocess.run(['bcftools', 'index', '-f', snp_file])
        if process.returncode != 0:
            print("VCF unsorted \nSorting..")
            subprocess.run(['bcftools', 'sort', snp_file, "-Oz", "-o", snp_file])
            subprocess.run(['bcftools', 'index', '-f', snp_file])       

    def set_chrom(snp_file, switch=True):
        """Description:
        sets variable to amount of chromosomes"""

        out = subprocess.run(['bcftools', 'query', '-f', '%CHROM\n', snp_file], stdout=subprocess.PIPE, text=True)
        ls = (out.stdout).split()
        ls_set = set(ls)
        if switch == True:
            try:
                print(f'Chromosomes: {listtostring(sorted(ls_set), ", ")}')
            except Exception:
                pass
        return len(ls_set), ls_set

    def check_chrom(snp_file, chr):

        out = subprocess.run(['bcftools', 'query', '-f', '%CHROM\n', snp_file], stdout=subprocess.PIPE, text=True)
        ls = (out.stdout).split()
        ls_set = set(ls)

        if chr != None:
            chr = list(dict.fromkeys(chr))
            chr_set = set(chr)
            chr_num = len(chr)

            if chr_set.issubset(ls_set) == True:
                pass
            else:
                sys.exit("Error: please provide the chromosome in the same format as in the VCF file")
        else:
            chr = sorted(ls_set)
            chr_num = len(ls_set)

        return chr, chr_num


    def filter_snps(min_af, subset, subset2, chr):
        """Description:
        Filters out SNPs with allele frequency below threshold via bcftools"""

        if chr == None:
            filtered = subprocess.run(['bcftools', 'view', '-m2', '-M2', '-v', 'snps', '-q', str(min_af), subset, '-Oz', '-o', subset2], stdout=subprocess.PIPE, text=True) 
            filtered.stdout
        else:
            filtered = subprocess.run(['bcftools', 'view', "-r", chr, '-m2', '-M2', '-v', 'snps', '-q', str(min_af), subset, '-Oz', '-o', subset2], stdout=subprocess.PIPE, text=True) 
            filtered.stdout

    def make_bed(subset2, chrom, memory, threads, list1):
        """Description:
        converts VCF file to PLINK BED files via plink"""

        string = "_"
        list2 = [l for l in list1 if string in l]
        if list2 != None:
            if chrom <= 24:
                make_bed = subprocess.run(['plink', '--vcf', f'{subset2}.vcf.gz', '--make-bed', '--out', subset2, '--mind', '1', '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--double-id', '--memory', str(memory),'--threads', str(threads)], stdout=subprocess.PIPE, text=True)
                make_bed.stdout            
            else:    
                make_bed = subprocess.run(['plink', '--vcf', f'{subset2}.vcf.gz', '--make-bed', '--out', subset2, '--mind', '1', '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--double-id', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)], stdout=subprocess.PIPE, text=True)
                make_bed.stdout
        else:
            #make_bed = subprocess.run(['plink', '--dummy', '15000', '2000000', '--make-bed', '--out', subset2, '--mind', '1', '--set-missing-var-ids', '@:#', '--allow-extra-chr', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)], stdout=subprocess.PIPE, text=True)
            if chrom <= 24:            
                make_bed = subprocess.run(['plink', '--vcf', f'{subset2}.vcf.gz', '--make-bed', '--out', subset2, '--mind', '1', '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--memory', str(memory),'--threads', str(threads)], stdout=subprocess.PIPE, text=True)            
                make_bed.stdout                            
            else:
                make_bed = subprocess.run(['plink', '--vcf', f'{subset2}.vcf.gz', '--make-bed', '--out', subset2, '--mind', '1', '--set-missing-var-ids', '@:#', '--allow-extra-chr', '0', '--memory', str(memory),'--threads', str(threads), '--chr-set', str(chrom)], stdout=subprocess.PIPE, text=True)
                make_bed.stdout

    def remove_files(subset, File, subset2, snp_file):
        """Description:
        removes no longer needed files"""

        try:
            os.remove(f'sub_{File}')
        except:
            pass
        try:
            os.remove(f'{subset}.txt')
        except:
            pass
        try:
            os.remove(f'{subset}.vcf.gz')
        except:
            pass
        try:
            os.remove(f'{subset2}.vcf.gz')
        except:
            pass

    def remove_covar(covar):
        """Description:
        removes temporary covariates file"""
        
        if covar != None:
            try:
                os.remove(f'sub_{covar}')
            except:
                pass


class Gemma:
    """Description:
    contains functions peforming GWAS analysis by calling GEMMA"""

    def write_returncodes(code):

        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"), 'a')
        file.write(str(code))
        file.close()
    
    def rel_matrix(prefix, Log, covar_file_name, model='-gk', n='1'):
        """Description:
        creates relatedness matrix and sets filename variable"""

        Log.print_log('Creating relatedness matrix..')
        # comment out following lines for testing if file is already in output
        if covar_file_name == None:
            process = subprocess.run(['gemma', '-bfile', prefix, model, n, '-o', prefix, '-outdir', "."])
        else:
            process = subprocess.run(['gemma', '-bfile', prefix, model, n, '-c', covar_file_name, '-o', prefix, '-outdir', "."])
        Gemma.write_returncodes(process.returncode)
        if process.returncode != 0:
            Log.print_log(f'Error: GEMMA was not able to complete the analysis')
            sys.exit(1)
        if n == '1':
            filename = f'{prefix}.cXX.txt'
        if n == '2':
            filename = None
        Log.print_log("Relatedness matrix created successfully")
        return filename

    def lm(prefix, prefix2, model, n, N, path, Log, covar_file_name):
        """Description:
        performs GWAS with linear model"""

        Log.print_log('Calculating linear model..')
        if covar_file_name == None:
            process = subprocess.run(['gemma', '-bfile', prefix, model, n, '-n', N, '-o', prefix2, '-outdir', path])
        else:
            process = subprocess.run(['gemma', '-bfile', prefix, model, n, '-n', N, '-c', covar_file_name, '-o', prefix2, '-outdir', path])
        Gemma.write_returncodes(process.returncode)
        if process.returncode != 0:
            Log.print_log(f'Error: GEMMA was not able to complete the analysis')
            sys.exit(1)
        Log.print_log("Linear model calculated successfully")

    def eigen(prefix, filename, model, Log):
        """Description:
        performs eigen-decomposition of relatedness matrix"""

        Log.print_log('Decomposing relatedness matrix..')
        process = subprocess.run(['gemma', '-bfile', prefix, '-k', filename, model, '-o', prefix, '-outdir', "."])
        Gemma.write_returncodes(process.returncode)
        if process.returncode != 0:
            Log.print_log(f'Error: GEMMA was not able to complete the analysis')
            sys.exit(1)
        Log.print_log("Eigen-decomposition of relatedness matrix successful")

    def lmm(pca, prefix, prefix2, filename, filename2, model, n, N, path, Log, covar_file_name):
        """Description:
        performs GWAS with linear mixed model"""

        Log.print_log('Calculating linear mixed model..')
        if covar_file_name == None:
            if pca != None:
                process = subprocess.run(['gemma', '-bfile', prefix, '-d', filename2, '-u', filename, model, n, '-n', N, '-o', prefix2, '-outdir', path], check=True)
            else:
                process = subprocess.run(['gemma', '-bfile', prefix, '-k', filename, model, n, '-n', N, '-o', prefix2, '-outdir', path], check=True)
        else:
            if pca != None:
                process = subprocess.run(['gemma', '-bfile', prefix, '-d', filename2, '-u', filename, model, n, '-n', N, "-c", covar_file_name, '-o', prefix2, '-outdir', path], check=True)
            else:
                process = subprocess.run(['gemma', '-bfile', prefix, '-k', filename, model, n, '-n', N, "-c", covar_file_name,'-o', prefix2, '-outdir', path], check=True)
        Gemma.write_returncodes(process.returncode)
        if process.returncode != 0:
            Log.print_log(f'Error code: \n{process.returncode} \nPossibly not enough memory available to process files.')
            sys.exit(1)
        Log.print_log("Linear mixed model calculated successfully")

    def bslmm(prefix, prefix2, model, n, N, path, Log, burn, sampling, snpmax):
        """Description:
        performs GWAS with bayesian sparse linear mixed model"""

        Log.print_log('Calculating bayesian sparse linear mixed model..')
        process = subprocess.run(['gemma', '-bfile', prefix, model, n, '-n', N, '-w', burn, '-s', sampling, '-smax', snpmax, '-o', prefix2, '-outdir', path], check=True)
        Gemma.write_returncodes(process.returncode)
        if process.returncode != 0:
            Log.print_log(f'Error: GEMMA was not able to complete the analysis')
            sys.exit(1)
        Log.print_log("Bayesian sparse linear mixed model calculated successfully")

    def run_gemma(prefix, prefix2, model, n, N, path, Log, filename, filename2, pca, covar_file_name, i, burn, sampling, snpmax):
        """Description:
        runs GEMMA dependent on input"""

        Log.print_log(f'Starting with {i} analysis..')
        if model not in ("-gk", "-eigen"):
            Log.print_log(f'Output will be saved in {path}/')

        if model == "-lm":
            Gemma.lm(prefix, prefix2, model, n, N, path, Log, covar_file_name)

        elif model == "-gk":
            Gemma.rel_matrix(prefix, Log, covar_file_name, model, n)

        elif model == "-eigen":
            if filename != None:
                Log.print_log("Reading relatedness matrix..")
            Gemma.eigen(prefix, filename, model, Log)

        elif model == "-lmm":
            if filename != None:
                if pca != None:
                    Log.print_log("Reading eigenvalues and eigenvectors..")
                else:
                    Log.print_log("Reading relatedness matrix..")
            Gemma.lmm(pca, prefix, prefix2, filename, filename2, model, n, N, path, Log, covar_file_name)

        elif model == "-bslmm":
            Gemma.bslmm(prefix, prefix2, model, n, N, path, Log, burn, sampling, snpmax)

        else:
            Log.print_log("No model was chosen!")

        Log.print_log(f'GEMMA executed successfully on {i}')


class Post_analysis:
    """Description:
    contains functions and subclasses regarding analysis of GWAS output"""

    def check_return_codes():

        code_file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_process_report.txt"), 'r')
        code_file_str = code_file.read()
        if "0" not in code_file_str:
            code_file.close()
            sys.exit(1)
        else:
            code_file.close()

    def load_df(prefix, case, path):
        """Description:
        reads GEMMA output file"""

        return pd.read_csv(os.path.join(path, f'{prefix}.{case}.txt'), sep='\t')

    def format_data(prefix, case, pcol, path, sep='\t', chromosome='chr'):
        """Description:
        formats dataframe for manhattan plot
        based on: https://github.com/Pudkip/Pyhattan/blob/master/Pyhattan/__init__.py"""

        data = pd.read_table(os.path.join(path, f'{prefix}.{case}.txt'), sep = sep)
        data = data[data[pcol] > 0.00000000000000000001]
        data['-log10(p_value)'] = -np.log10(data[pcol])
        data = data[data['-log10(p_value)'].notna()]
        data[chromosome] = data[chromosome].astype('category')
        data['ind'] = range(len(data))
        data_grouped = data.groupby((chromosome))
        return data, data_grouped

    def manh_plot(df, Log, prefix, pcol, path, sigval, x, nolabel, colors = ['#E24E42', '#008F95'], refSNP = False):
        """Description:
        creates manhattan plot from prepared dataframe and saves plot
        based on: https://github.com/Pudkip/Pyhattan/blob/master/Pyhattan/__init__.py"""

        Log.print_log(f'Creating Manhattan plot of {pcol}..')

        data = df[0]
        data_grouped = df[1]
        n = 0
        if data.empty == True:
            Log.print_log("GEMMA couldn't calculate any meaningful values!")
        else:
            timer = time.perf_counter()
            fig = plt.figure(figsize=(16,12))
            ax = fig.add_subplot(111)

            x_labels = []
            x_labels_pos = []
            texts = []

            for num, (name, group) in enumerate(data_grouped):
                group.plot(kind='scatter', x='ind', y='-log10(p_value)', color=colors[num % len(colors)], ax=ax, s=fontsize3)
                x_labels.append(name)
                x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels)
            ax.set_xlim([0, len(data)])
            ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
            ax.set_xlabel('Chromosome', fontsize=fontsize, fontweight="bold")
            ax.set_title("Manhattan plot", fontsize=fontsize)
            ax.set_ylabel(f'-log10({pcol})', fontsize=fontsize, fontweight="bold")
            plt.xticks(fontsize=fontsize2, rotation=45, fontweight="bold")
            plt.yticks(fontsize=fontsize2, fontweight="bold")
            np.random.seed(0)

            # Bonferroni correction
            sig_level = 0
            if sigval != None:
                sig_level = sigval
            if sigval == None and x != 0:
                sig_level = (0.05 / x)
                sigval = -np.log10(sig_level)

            if refSNP:
                for index, row in data.iterrows():
                    if sigval > 0:
                        if row['-log10(p_value)'] >= sigval:
                            n += 1
                            if nolabel == False:                  
                                texts.append(plt.text(index, row['-log10(p_value)']+np.random.random()/100, str(row[refSNP]), fontsize="medium", fontweight="bold"))
                if sigval > 0:
                    plt.axhline(y=sigval, color='black', linestyle='-', linewidth = fontsize6)
                    Log.print_log(f'Number of significant SNPs: {n} \nLevel of significance: {np.format_float_scientific(sig_level, precision=2)} \nNumber of total SNPs: {x}')
                    if nolabel == False:
                        adjust_text(texts, data.index.values, data['-log10(p_value)'].values, autoalign='y', ha='left', only_move={'text':'y'}, expand_text=(1.02, 1.02), expand_align=(1.02, 1.02), force_text=(0,0.7), lim=250, arrowprops=dict(arrowstyle="-", color='k', lw=0.5, alpha=0.6))

            file_path = os.path.join(path, "manhattan")
            make_dir(file_path)
            plt.savefig(os.path.join(file_path, f'{pcol}_manh_{prefix}.png'))
            plt.close()
            timer_end = time.perf_counter()
            timer_total = round(timer_end - timer, 2)
            Log.print_log(f'Manhattan plot saved as "{pcol}_manh_{prefix}.png" in {file_path} (Duration: {runtime_format(timer_total)})')

        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_snpcount_total.txt"), 'a')
        file.write(f'{x}\n')
        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_snpcount_sig.txt"), 'a')
        file.write(f'{n}\n')
        file = open(os.path.join("_vcf2gwas_temp", "vcf2gwas_sig_level.txt"), 'a')
        file.write(f'{np.format_float_scientific(sig_level, precision=2)}\n')
        file.close()

        return n

    def ppoints(n):
        """Description:
        generates sequence of probability points
        based on: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/ppoints"""

        a = 3/8 if n <= 10 else 1/2
        try:
            n = (len(n))
        except Exception as e:
            n = n
        l = []
        for i in range(n):
            l.append((i - a)/(n + 1 - 2*a))
        return l

    def qq_plot(df, pcol, prefix, path, Log):
        """Description:
        creates qq-plot from dataframe and saves plot
        based on: https://github.com/stephenturner/qqman/blob/master/R/qq.R"""

        Log.print_log('Creating QQ-plot..')
        df = df[df[pcol] > 0.0000000001]
        df['-log10(p_value)'] = -np.log10(df[pcol])
        df = df[df['-log10(p_value)'].notna()]
        if df.empty == True:
            Log.print_log("GEMMA couldn't calculate any meaningful values!")
        else:
            timer = time.perf_counter()
            a = df['-log10(p_value)'].values.tolist()
            a = sorted(a, reverse=True)
            a_max = (max(a)+0.5)

            c = len(a)

            b = Lin_models.ppoints(c)
            b = [-np.log10(i) if i > 0 else 0 for i in b]
            b = sorted(b, reverse=True)
            b_max = (np.nanmax(b)+0.5)
            q = np.maximum(a_max, b_max)

            fig, ax = plt.subplots(figsize=(16,12))
            ax.plot(b,a, ls="", marker="o", alpha=0.5, markeredgecolor="black", markeredgewidth=1)
            ax.set(xlim=(0, (b_max)), ylim=(0, a_max))
            x = np.linspace(0, q)
            ax.plot(x,x, color="k", ls="--")
            ax.set_title("QQ-Plot", fontsize=fontsize)
            ax.set_xlabel("Expected  "r'$-log_{10}(p)$', fontsize=fontsize, fontweight="bold")
            ax.set_ylabel("Observed  "r'$-log_{10}(p)$', fontsize=fontsize, fontweight="bold")
            plt.xticks(fontsize=fontsize2, fontweight="bold")
            plt.yticks(fontsize=fontsize2, fontweight="bold")

            file_path = os.path.join(path, "QQ")
            make_dir(file_path)
            plt.savefig(os.path.join(file_path, f'{pcol}_qq_{prefix}.png'))
            plt.close()
            timer_end = time.perf_counter()
            timer_total = round(timer_end - timer, 2)
            Log.print_log(f'QQ-plot saved as "{pcol}_qq_{prefix}.png" in {file_path} (Duration: {runtime_format(timer_total)})')

    def make_top_list(df, top_list, n, x):
        """Description:
        returns list of the top n SNPs with highest p-value"""

        if x > n:
            n = x
        new_df = df.head(n)
        new_df2 = new_df.dropna(axis=0, how="any")
        for i in ["chr", "rs", "ps"]:
            col = new_df[i].astype(str).tolist()
            top_list.append(col)

    def print_top_list(l, cols, path, pc_prefix, snp_prefix):
        """Description:
        saves dataframe with top SNPs as file"""

        df = pd.DataFrame(l)
        df = df.T
        iterables = [cols, ["chr", "SNP_ID", "SNP_pos"]]
        multcols = pd.MultiIndex.from_product(iterables)
        df.columns = multcols
        df.to_csv(os.path.join(path, f'top_SNPs{pc_prefix}_{snp_prefix}.csv'), sep=',')

    def run_postprocessing(top_ten, Log, model, n, prefix2, path, n_top, i, sigval, nolabel):
        """Description:
        runs post-processing dependent on input"""

        Log.print_log(f'Starting post-processing of {i}..')
        if model in ("-lm", "-lmm"):

            pcol = Lin_models.set_pcol(n)
            
            Log.print_log("Summarizing p-values..")
            for p in pcol:
                df = Lin_models.load_df(prefix2, "assoc", path)
                x = len(df.index)
                df2 = Lin_models.get_p_values(df, p, prefix2, path)
                Log.print_log(f'Variants with the best {p} score saved in {os.path.join(path, "best_p-values")}')

                df3 = Lin_models.format_data(prefix2, "assoc", p, path)
                n = Lin_models.manh_plot(df3, Log, prefix2, p, path, sigval, x, nolabel, refSNP="rs")
                Lin_models.qq_plot(df, p, prefix2, path, Log)

            Lin_models.make_top_list(df2, top_ten, n_top, n)

        elif model == "-bslmm":
            # procedure based on steps outlined in: 
            # http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf

            ## get hyperparameters
            Log.print_log("Summarizing hyperparameters..")
            # format column names
            Bslmm.format_col(prefix2, path)
            # Load hyperparameters
            df = Bslmm.load_df(prefix2, "hyp", path)
            df = Bslmm.rm_unnamed(df)
            # Get mean, median, and 95% ETPI of hyperparameters
            a = Bslmm.get_hyper(df, prefix2, path)
            Log.print_log(f'Mean, median, and 95% ETPI of hyperparameters saved as "hyperparameters.csv" in {os.path.join(path, "hyperparameters")}')
            # plot traces and distributions of hyperparameters
            Bslmm.diagnostics(df, a, prefix2, path, Log)

            ## get parameters
            Log.print_log("Summarizing parameters..")
            # load parameters
            df2 = Bslmm.load_df(prefix2, "param", path)
            x = len(df.index)
            # Get variants with sparse effect size on phenotypes
            Bslmm.get_eff(df2, prefix2, path)
            #df3 = Bslmm.get_eff(df2, prefix2, path)
            Log.print_log(f'Variants with the highest sparse effects saved in {os.path.join(path, "highest_effects")}')
            # Get variants with high Posterior Inclusion Probability (PIP) == gamma
            #Bslmm.get_pip(df2, prefix2, path)
            df3 = Bslmm.get_pip(df2, prefix2, path)
            Log.print_log(f'Variants with high Posterior Inclusion Probability (PIP) (== gamma) saved in {os.path.join(path, "high_PIP")}')

            # plot variants PIPs across linkage groups/chromosomes
            df4 = Bslmm.format_data(prefix2, "param", "gamma", path)
            n1 = Bslmm.manh_plot(df4, Log, prefix2, "gamma", path, sigval, x, nolabel, refSNP="rs")
            # plot effect sizes
            df5 = Bslmm.format_data(prefix2, "param", "eff", path)
            n2 = Bslmm.manh_plot(df5, Log, prefix2, "eff", path, sigval, x, nolabel, refSNP="rs")
            n = max([n1, n2])

            Bslmm.make_top_list(df3, top_ten, n_top, n)
        
        else:
            Log.print_log("No post-processing necessary!")
        Log.print_log(f'Analysis of {i} completed\n')

class Lin_models(Post_analysis):
    """Description:
    Post-analysis subclass containing functions specific to analysis with linear models"""

    def set_pcol(n):
        """Description:
        sets variable to correct p-value column"""

        pdict = {"1": ["p_wald"], "2": ["p_lrt"], "3": ["p_score"], "4": ["p_lrt", "p_score", "p_wald"]}
        return pdict[str(n)]

    def get_p_values(df, pcol, prefix, path):
        """Description:
        sorts dataframe by p-value and save files with top variants"""

        # variants with best p-values
        df2 = df.sort_values(by=pcol, ascending=True)
        df3 = df
        df3 = df3[df3[pcol].notna()]

        # top 1% variants (below 1% quantile)
        top1 = df2[df2[pcol]<df2[pcol].quantile(0.01)]
        # top 0.1% variants (below 0.1% quantile)
        top01 = df2[df2[pcol]<df2[pcol].quantile(0.001)]
        # top 0.01% variants (below 0.01% quantile)
        top001 = df2[df2[pcol]<df2[pcol].quantile(0.0001)]
        #write files
        file_path = os.path.join(path, "best_p-values")
        make_dir(file_path)
        top1.to_csv(os.path.join(file_path, f'{pcol}_{prefix}_top1.csv'))
        top01.to_csv(os.path.join(file_path, f'{pcol}_{prefix}_top01.csv'))
        top001.to_csv(os.path.join(file_path, f'{pcol}_{prefix}_top001.csv'))
        if df3.empty == True:
            df2["rs"] = "NaN"
            df2["chr"] = "NaN"
            df2["ps"] = "NaN"
        return df2

class Bslmm(Post_analysis):
    """Description:
    Post-analysis subclass containing functions specific to analysis with bayesian model""" 

    def format_col(prefix, path):
        """Description:
        formats column names in GEMMA output file"""

        with open(os.path.join(path, f'{prefix}.hyp.txt')) as f:
            newText=f.read().replace('h ', 'h').replace(' pve', 'pve').replace(' rho', 'rho').replace(' pge', 'pge').replace(' pi', 'pi').replace(' n_gamma', 'n_gamma\t')
        with open(os.path.join(path, f'{prefix}.hyp.txt'), "w") as f:
            f.write(newText)

    def rm_unnamed(df):
        """Description:
        removes column"""

        return df.loc[:, ~df.columns.str.contains('^Unnamed')]

    def get_hyper(df, prefix, path):
        """Description:
        prints hyperparameters to file"""

        a = df.columns
        b = df.mean(axis=0)
        df3_temp = pd.DataFrame(list(zip(a,b)))
        c = df.quantile([0.5, 0.025,0.975]).T
        c = c.reset_index(drop=True)
        df3 = pd.concat([df3_temp,c], axis=1)
        df3.columns = ["hyperparam", "mean", "median", "2.5%", "97.5"]
        # save hyperparameters as csv
        make_dir(os.path.join(path, "hyperparameters"))
        df3.to_csv(os.path.join(path, "hyperparameters", f'{prefix}_hyperparameters.csv'))
        return a

    def diagnostics(df, a, prefix, path, Log):
        """Description:
        saves diagnostic plots of hyperparameters"""

        Log.print_log("Plotting traces and distributions of hyperparameters..")
        make_dir(os.path.join(path, "diagnostics"))
        for i in a:
            #create grid
            fig =plt.figure(figsize=(16,12))
            grid = plt.GridSpec(2,2, wspace=0.4, hspace=0.3)
            line = fig.add_subplot(grid[0,:2])
            plt.xlabel("Index", fontsize=fontsize, fontweight="bold")
            plt.ylabel(i, fontsize=fontsize, fontweight="bold")
            plt.xticks(fontsize=fontsize3, fontweight="bold")
            plt.yticks(fontsize=fontsize3, fontweight="bold")
            hist = fig.add_subplot(grid[1,0])
            plt.xlabel(i, fontsize=fontsize, fontweight="bold")
            plt.ylabel("Frequency", fontsize=fontsize, fontweight="bold")
            plt.xticks(fontsize=fontsize3, fontweight="bold")
            plt.yticks(fontsize=fontsize3, fontweight="bold")
            dens = fig.add_subplot(grid[1,1])
            plt.xlabel(i, fontsize=fontsize, fontweight="bold")
            plt.ylabel("Density", fontsize=fontsize, fontweight="bold")
            plt.xticks(fontsize=fontsize3, fontweight="bold")
            plt.yticks(fontsize=fontsize3, fontweight="bold")
            fig.suptitle(f'Diagnostic plots of {i}', fontsize=fontsize)
            #add datapoints
            sns.lineplot(data=df[i], ax=line)
            sns.histplot(df[i], ax=hist, kde=False)
            try:
                sns.kdeplot(df[i], ax=dens)
            except Exception as e:
                pass
            #save figures  
            plt.savefig(os.path.join(path, "diagnostics", f'{prefix}_{i}_hyperparameter.png'))
            plt.close()
            Log.print_log(f'Plot of {i} saved in {os.path.join(path, "diagnostics")}')

    def get_eff(df, prefix, path):
        """Description:
        calculates sparse effect size, sorts and prints to file"""

        # add sparse effect size (= beta * gamma) to data frame
        eff = df.beta * df.gamma
        df["eff"] = eff

        # sort by decreasing effect size
        df5 = df.sort_values(by="eff", ascending=False)
        # get variants with effect size > 0
        df6 = df5[df5["eff"] > 0]
        # variants with the highest sparse effects 
        # top 1% variants (above 99% quantile)
        top1 = df6[df6["eff"]>df6["eff"].quantile(0.99)]
        # top 0.1% variants (above 99.9% quantile)
        top01 = df6[df6["eff"]>df6["eff"].quantile(0.999)]
        # top 0.01% variants (above 99.99% quantile)
        top001 = df6[df6["eff"]>df6["eff"].quantile(0.9999)]
        #write files
        df.to_csv(os.path.join(path, f'{prefix}.param.txt'), sep='\t', index=False)
        file_path = os.path.join(path, "highest_effects")
        make_dir(file_path)
        top1.to_csv(os.path.join(file_path, f'{prefix}_top1eff.csv'))
        top01.to_csv(os.path.join(file_path, f'{prefix}_top01eff.csv'))
        top001.to_csv(os.path.join(file_path, f'{prefix}_top001eff.csv'))
        if df6.empty == True:
            df5["rs"] = "NaN"
        return df5

    def get_pip(df, prefix, path):
        """Description:
        sorts variants by descending posterior inclusion propability and prints to file"""

        # sort variants by descending PIP
        df7 = df.sort_values(by="gamma", ascending=False)
        df8 = df
        df8 = df8[df8["gamma"].notna()]
        # sets of variants above a certain threshold
        # variants with effect in 1% MCMC samples or more
        pip01 = df7[df7["gamma"]>0.01]
        # variants with effect in 10% MCMC samples or more
        pip10 = df7[df7["gamma"]>0.1]
        # variants with effect in 25% MCMC samples or more
        pip25 = df7[df7["gamma"]>0.25]
        # variants with effect in 50% MCMC samples or more
        pip50 = df7[df7["gamma"]>0.5]
        #write files
        file_path = os.path.join(path, "high_PIP")
        make_dir(file_path)
        pip01.to_csv(os.path.join(file_path, f'{prefix}_pip01.csv'))
        pip10.to_csv(os.path.join(file_path, f'{prefix}_pip10.csv'))
        pip25.to_csv(os.path.join(file_path, f'{prefix}_pip25.csv'))
        pip50.to_csv(os.path.join(file_path, f'{prefix}_pip50.csv'))
        if df8.empty == True:
            df7["rs"] = "NaN"
        return df7


class Summary:
    """Description:
    contains functions regarding summary of GWAS output"""

    def summarizer(path, path2, pc_prefix, snp_prefix, n, Log, prefix_list):
        """Description:
        reads all top_SNPs files and summarizes them"""

        Log.print_log("Summarizing top SNPs..")
        # read files as dataframes and set as variables
        filenames = []
        dfnames = []
        x = 0
        for prefix in prefix_list:
            for file in os.listdir(path):
                if file == f'top_SNPs{prefix}_{snp_prefix}.csv':
                    x = x + 1 
                    filename = os.path.join(path, file)
                    globals()[f'df{x}'] = pd.read_csv(filename, header=[0,1], index_col=0, sep=',')
                    temp = list(eval(f'df{x}').columns.get_level_values(1))
                    globals()[f'df{x}'].loc[(n)] = temp
                    dfnames.append(f'df{x}')
        if dfnames == []:
            Log.print_log("Couldn't find files to summarize!")
            values = None
        else:
            # concat the dataframes
            if len(dfnames) == 1:
                dfx = df1
            else:
                for i in range(len(dfnames)-1):
                    a = i+1
                    b = i+2
                    globals()[f'df{b}'] = pd.concat([eval(f'df{a}'),eval(f'df{b}')])
                    dfx = eval(f'df{b}')
            # make lists from snp, pos and chr columns and remove nan
            cols = dfx.columns
            num = int(len(cols))
            x1 = np.arange(0, num, 3)
            x2 = np.arange(1, (num + 1), 3)
            x3 = np.arange(2, (num + 2), 3)
            l2 = []
            l3 = []
            l4 = []
            for i in x1:
                temp = dfx.iloc[:,i].astype(str).tolist()
                l2.extend(temp)
            for i in x2:
                temp = dfx.iloc[:,i].astype(str).tolist()
                l3.extend(temp)
            for i in x3:
                temp = dfx.iloc[:,i].astype(str).tolist()
                l4.extend(temp)
            l2 = [i for i in l2 if i != "nan"]
            l3 = [i for i in l3 if i != "nan"]
            l4 = [i for i in l4 if i != "nan"]
            for i in [l2,l3,l4]:
                if "chr" in i:
                    l2_corr = i
                    l2_corr = [x for x in i if x != "chr"]
                if "SNP_ID" in i:
                    l3_corr = i
                    l3_corr = [x for x in i if x != "SNP_ID"]
                if "SNP_pos" in i:
                    l4_corr = i
                    l4_corr = [x for x in i if x != "SNP_pos"]
            # concat the lists into an array, convert to dataframe, set types and columns names
            l = np.array([l3_corr,l4_corr,l2_corr], dtype=object)
            values = pd.DataFrame(l.T)
            values[1] = values[1].apply(pd.to_numeric)
            values[1] = values[1].astype(int)
            values.columns = ["SNP_ID", "SNP_pos", "chr"]
            if values.empty == True:
                    Log.print_log(str("No SNPs in the summarized files!"))
            values_list = []
            # remove unwanted SNPs
            if len(x1) == 1 and len(dfnames) == 1:
                values = values.dropna()
                values_list.append(values)
            else:
                # count occurrences of SNPs and filter
                values = values.value_counts()
                values = pd.DataFrame(values)
                values = values.reset_index()
                values.columns = ["SNP_ID", "SNP_pos", "chr", "count"]
                # remove unwanted SNPs
                values2 = values.dropna()
                values_list.append(values2)
                values = values.where(values["count"]>1)
                values = values.dropna()
                values[["SNP_pos", "count"]] = values[["SNP_pos", "count"]].astype(int)       
                try:
                    # plot SNP counts
                    y = values["count"]
                    fig = plt.figure(figsize=(16,12))
                    sns.scatterplot(data=values, x="SNP_ID", y="count", s=(fontsize*2))
                    sns.despine(offset=10)
                    plt.xticks(rotation=45, fontsize=fontsize4, fontweight="bold")
                    plt.xlabel("SNP ID", fontsize=fontsize, fontweight="bold")
                    plt.yticks(np.arange(0, max(y)+2, 1), fontsize=fontsize2, fontweight="bold")
                    plt.ylabel("Occurrence of SNP", fontsize=fontsize, fontweight="bold")
                    plt.savefig(os.path.join(path2, f'summarized_top_SNPs{pc_prefix}_{snp_prefix}.png'))
                    plt.close()
                    values_list.append(values)
                except Exception as e:
                    Log.print_log("No SNP occurred more than once")
            # make list of phenoytpes where SNPs occurred and add as column
            addstring = ["_complete", ""]
            c = 0
            for values in values_list:
                names = []
                for i in values.index:
                    xy = values.iloc[i,0]
                    dfx2 = dfx.where(dfx.isin([xy]))
                    name = dfx2.dropna(axis=1, how="all")
                    name = name.columns.get_level_values(0).astype(str).to_list()
                    name = listtostring(name, ', ')
                    names.append(name)
                values["phenotypes"] = names
                #save as file
                if len(dfnames) == 1:
                    filename = f'summarized_top_SNPs{addstring[c]}{pc_prefix}_{snp_prefix}.csv'
                    values.to_csv(os.path.join(path2, filename) , index=False, sep=',')
                    Log.print_log(f'Top SNPs summarized and saved as "{filename}" in {path2}')
                else:
                    filename = f'summarized_top_SNPs{addstring[c]}_{snp_prefix}.csv'
                    values.to_csv(os.path.join(path2, filename), index=False, sep=',')
                    Log.print_log(f'Top SNPs summarized and saved as "{filename}" in {path2}')
                c += 1
                filenames.append(os.path.join(path2, filename))
        return filenames

    def gff_converter(f):

        df = pd.read_csv(f, sep="\t", header=0, names=["chr", "source", "type", "start", "stop", "score", "strand", "phase", "attributes"], comment="#")
        df.where(df["type"]=="gene", inplace=True)
        df.dropna(how="any", inplace=True)
        df2 = pd.DataFrame(df.attributes.str.split(";").tolist())
        l_ID = []
        l_note = []
        l_name = []

        for string, l in zip(["ID", "Note", "Name"], [l_ID, l_note, l_name]):
            for row in range(len(df2.index)):
                val = "NaN"
                for col in df2.columns:
                    try:
                        if df2.iloc[row,col].startswith(string):
                            val = df2.iloc[row,col]
                            val = val.removeprefix(f'{string}=')
                    except:
                        pass
                l.append(val)
        
        l_note2 = [i for i in l_note if str(i) != "NaN"]
        if l_note2 == []:
            l_descr = []
            for row in range(len(df2.index)):
                val = "NaN"
                for col in df2.columns:
                    try:
                        if df2.iloc[row,col].lower().startswith("description"):
                            val = df2.iloc[row,col]
                            try:
                                val = val.removeprefix(f'description=')
                            except:
                                val = val.removeprefix(f'Description=')
                    except:
                        pass
                l_descr.append(val)
            l_note = l_descr

        df["ID"] = l_ID
        df["comment"] = l_note
        df["name"] = l_name
        df.drop("attributes", axis=1, inplace=True)
        df.replace("NaN", np.nan, inplace=True)
        df[["start", "stop"]] = df[["start", "stop"]].round().astype("int")
        df = df.reset_index(drop=True)
        return df

    def chr_converter(df, chr):

        df["chr"] = df["chr"].astype("str")
        chr2 = df["chr"].tolist()
        chr_dict = {}
        chr_set = set(chr)
        chr2_set = set(chr2)

        if len(chr_set.intersection(chr2_set)) != 0:
            pass
        else:
            keys = []
            values = []
            for elem in chr_set:
                for c in chr2_set:
                    if elem.lower() in c.lower():
                        keys.append(c)
                        values.append(elem)
                    elif c.lower() in elem.lower():
                        keys.append(c)
                        values.append(elem)
            for i in range(len(keys)):
                chr_dict[keys[i]] = values[i]
            if not chr_dict:
                sys.exit("Error: Please make sure that the chromosomes in the VCF and gene file have the same format")

        df["chr"].replace(chr_dict, inplace=True)
        df = df[df["chr"].isin(chr_set)]
        return df

    def gene_compare(filenames, gene_file, gene_file_path, gene_thresh, path, pc_prefix, snp_prefix, chr_list, Log):
        """Description:
        takes dataframe of summarized SNPs, calculates distances to genes in gene input file and saves results in new file.
        Columns in gene file: start, stop, chr, name (optional), ID (optional), comment (optional)"""
        
        Log.print_log("Comparing top SNPs to genes..")
        addstring = ["_complete", ""]
        c = 0        
        for filename in filenames:
            values = pd.read_csv(filename)
            if values.empty == True:
                Log.print_log("No SNPs present to compare to genes")
            else:
                n = 0
                for file in os.listdir(path):
                    if file == f'summarized_top_SNPs{pc_prefix}_{snp_prefix}.csv':
                        n += 1
                #read gene .csv file
                if gene_file_path.endswith(".gff"):
                    df = Summary.gff_converter(gene_file_path)
                    gene_file_name = gene_file.removesuffix(".gff")
                else:
                    df = pd.read_csv(gene_file_path, sep=',')
                    gene_file_name = gene_file.removesuffix(".csv")
                df = Summary.chr_converter(df, chr_list)
                #set list variables
                gene_dist_left = []
                gene_dist_right = []
                gene_ID_left = []
                gene_ID_right = []
                gene_annot_left = []
                gene_annot_right = []
                gene_comment_left = []
                gene_comment_right = []
                #check file
                try:
                    df_new = df[["chr", "start", "stop"]]
                    df_new = df_new.dropna(subset=["start", "stop"])
                    df_new[["start", "stop"]] = df_new[["start", "stop"]].astype("int32")
                    check = True
                except Exception:
                    check = False
                if check == True:
                    #loop through columns of gene file and calculate distances between SNPs and genes
                    for x in values.index:
                        df_new = df[["chr", "start", "stop"]]
                        try:
                            df_new["ID"] = df["ID"]
                        except Exception:
                            pass
                        try:
                            df_new["name"] = df["name"]
                        except Exception:
                            pass
                        try:
                            df_new["comment"] = df["comment"]
                        except Exception:
                            pass
                        df_new["Chr"] = values.loc[x,"chr"]
                        df_new[["chr", "Chr"]] = df_new[["chr", "Chr"]].astype("str")
                        df_new["pos"] = values.loc[x,"SNP_pos"]
                        df_new["SNP_ID"] = values.loc[x,"SNP_ID"]
                        df_new = df_new.dropna(subset=["start", "stop"])
                        df_new[["pos", "start", "stop"]] = df_new[["pos", "start", "stop"]].astype("int32")
                        df_new = df_new[df_new.chr==df_new.Chr]
                        df_new["POSleft"] = df_new["pos"] - df_new["stop"]
                        df_new["POSright"] = df_new["start"] - df_new["pos"]
                        df_new[["POSleft", "POSright"]] = df_new[["POSleft", "POSright"]].astype("int32")
                        #filter results and append to lists
                        min_value1 = df_new[(df_new["POSleft"] > 0) & (df_new["POSleft"] < gene_thresh)]
                        min_value1 = min_value1[min_value1.POSleft == min_value1.POSleft.min()]
                        try:
                            gene_dist_left.append(min_value1.iloc[0]["POSleft"])
                        except Exception:
                            gene_dist_left.append(np.nan)
                        try:
                            gene_ID_left.append(min_value1.iloc[0]["ID"])
                        except Exception:
                            gene_ID_left.append(np.nan)
                        try:
                            gene_annot_left.append(min_value1.iloc[0]["name"])
                        except Exception:
                            gene_annot_left.append(np.nan)
                        try:
                            gene_comment_left.append(min_value1.iloc[0]["comment"])
                        except Exception:
                            gene_comment_left.append(np.nan)
                        min_value2 = df_new[(df_new["POSright"] > 0) & (df_new["POSright"] < gene_thresh)]
                        min_value2 = min_value2[min_value2.POSright == min_value2.POSright.min()]
                        try:
                            gene_dist_right.append(min_value2.iloc[0]["POSright"])
                        except Exception:
                            gene_dist_right.append(np.nan)
                        try:
                            gene_ID_right.append(min_value2.iloc[0]["ID"])
                        except Exception:
                            gene_ID_right.append(np.nan)
                        try:
                            gene_annot_right.append(min_value2.iloc[0]["name"])
                        except Exception:
                            gene_annot_right.append(np.nan)
                        try:
                            gene_comment_right.append(min_value2.iloc[0]["comment"])
                        except Exception:
                            gene_comment_right.append(np.nan)
                    # add lists to dataframe
                    pos_col = values["SNP_pos"].to_list()
                    values = values.drop("SNP_pos", axis=1)
                    values["gene_ID(up)"] = gene_ID_left
                    values["gene_comment(up)"] = gene_comment_left
                    values["gene_name(up)"] = gene_annot_left
                    values["gene_distance(up)"] = gene_dist_left
                    values["SNP_pos"] = pos_col
                    values["gene_distance(down)"] = gene_dist_right
                    values["gene_name(down)"] = gene_annot_right
                    values["gene_comment(down)"] = gene_comment_right
                    values["gene_ID(down)"] = gene_ID_right
                    # sort and remove empty columns
                    values = values.sort_values(by=["gene_distance(up)","gene_distance(down)"])
                    values.dropna(subset=["gene_distance(up)", "gene_distance(down)"], how="all", inplace=True)
                    values.dropna(axis=1, how="all", inplace=True)
                    try:
                        values["gene_distance(up)"] = values["gene_distance(up)"].astype("Int32")
                    except Exception:
                        pass
                    try:
                        values["gene_distance(down)"] = values["gene_distance(down)"].astype("Int32")
                    except Exception:
                        pass
                    try:
                        values["count"] = values["count"].astype("int32")
                    except Exception:
                        pass
                    #make multicolumns
                    if len(values.columns) == 12:
                        multicols = [
                            np.array(["", "", "", "Upstream gene", "Upstream gene", "Upstream gene", "Upstream gene", "", "Downstream gene", "Downstream gene", "Downstream gene", "Downstream gene"]),
                            np.array(["SNP ID", "Chr", "Phenotype", "ID", "Comment", "Name", "Distance", "SNP position", "Distance", "Name", "Comment", "ID"])
                        ]
                        values.columns = multicols
                    elif len(values.columns) == 13:
                        multicols = [
                            np.array(["", "", "", "", "Upstream gene", "Upstream gene", "Upstream gene", "Upstream gene", "", "Downstream gene", "Downstream gene", "Downstream gene", "Downstream gene"]),
                            np.array(["SNP ID", "Chr", "Count", "Phenotype", "ID", "Comment", "Name", "Distance", "SNP position", "Distance", "Name", "Comment", "ID"])
                        ]
                        values.columns = multicols                       
                    #save file
                    if values.empty:
                        Log.print_log("Info: No SNPs could be compared to genes, please check the species or gene file selection")
                    else:
                        if n > 0:
                            values.to_csv(os.path.join(path, f'{gene_file_name}_compared_summarized_top_SNPs{addstring[c]}{pc_prefix}_{snp_prefix}.csv'), index=False)
                            Log.print_log(f'Top SNPs compared to genes and saved as "{gene_file_name}_compared_summarized_top_SNPs{addstring[c]}{pc_prefix}_{snp_prefix}.csv" in {path}')
                        else:
                            values.to_csv(os.path.join(path, f'compared_summarized_top_SNPs{addstring[c]}_{snp_prefix}.csv'), index=False)
                            Log.print_log(f'Top SNPs compared to genes and saved as "compared_summarized_top_SNPs{addstring[c]}_{snp_prefix}.csv" in {path}')
                else:
                    Log.print_log(f'Could not parse contents of {gene_file}.\nPlease provide the file in the right format.')
            c += 1