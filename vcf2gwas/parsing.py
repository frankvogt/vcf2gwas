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

import argparse
import os
import subprocess
import multiprocessing as mp

try:
    from psutil import virtual_memory
except ModuleNotFoundError:
    subprocess.run(["conda", "install", "-c", "anaconda", "psutil==5.8*"])
    from psutil import virtual_memory


def getArgs(argv=None):
    """Description:
    Sets up Argument Parser and returns input arguments"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='Command-line interface for vcf2gwas.\n \nExample usage: vcf2gwas -v <VCF file> -pf <phenotype file> -ap -lmm', epilog="For a detailed description of all options, please refer to the manual.")

    parser.add_argument('--version', action='version', version='%(prog)s 0.6.2')
    parser.add_argument("-v", "--vcf", metavar="<filename>", required=True, type=str, help="(required) Genotype .vcf or .vcf.gz filename")
    parser.add_argument("-pf", "--pfile", metavar="<filename>", action="append", type=str, help="specify phenotype filename \ncomma separated .csv file \nfirst column individuals, second column and onwards phenotypes")
    parser.add_argument("-cf", "--cfile", metavar="<filename>", type=str, help="specify covariate filename \ncomma separated .csv file \nfirst column individuals, second column and onwards covariates\n")
    parser.add_argument("-p", "--pheno", metavar="<int>", action="append", type=int, help="specify phenotypes used for analysis: \n '1' selects first phenotype from phenotype file (second column),\n '2' the second phenotype (third column) and so on")
    parser.add_argument("-c", "--covar", metavar="<int>", action="append", type=int, help="specify covariates used for analysis: \n '1' selects first covariate from covariate file (second column),\n '2' the second covariate (third column) and so on")
    parser.add_argument("-gf", "--genefile", metavar="<filename>", type=str, help= "specify gene file name \nresulting SNPs from analysis will be compared to these genes \ncomma separated .csv file \nfile must contain at least three columns: \n -'chr' column containing chromosome value (same format as in VCF file \n -'start' column containing start position of gene \n -'stop' column containing stop position of gene \nfor further formatting information, please refer to the manual")
    parser.add_argument("-gt", "--genethresh", metavar="<int>", type=int, default=100000, help="set gene distance threshold when comparing genes (in bp) (default: %(default)s) \nonly SNPs with distances below threshold will be considered for comparison of each gene")
    parser.add_argument("-q", "--minaf", metavar="<float>", type=float, default=0.01, help="minimum allele frequency of sites to be used (default: %(default)s) \ninput value needs to be a value between 0.0 and 1.0 ")
    parser.add_argument("-t", "--topsnp", metavar="<int>", type=int, default=15, help="number of top SNPs of each phenotype to be summarized (default: %(default)s) \nafter analysis the specified amount of top SNPs from each phenotype will be considered")
    parser.add_argument("-M", "--memory", metavar="<int>", type=int, default=int(((virtual_memory().total/1e9)//2)*1e3), help="set memory usage (in MB) \nif not specified, half of total memory will be used (%(default)s MB)")
    parser.add_argument("-T", "--threads", metavar="<int>", type=int, default=mp.cpu_count()-1, help="set core usage \nif not specified, all available logical cores minus 1 will be used (%(default)s cores)")
    parser.add_argument("-k", "--relmatrix", metavar="<filename>", type=str, help="specify relatedness matrix file name \nfor formatting details please refer to the manual")
    parser.add_argument("-P", "--PCA", metavar="<int>", type=int, nargs="?", const=2, help="perform PCA on phenotypes and use resulting PCs as phenotypes for GEMMA analysis \noptional: set amount of PCs to be calculated (default: %(const)s) \nrecommended amount of PCs: 2 - 10")
    parser.add_argument("-U", "--UMAP", metavar="<int>", type=int, nargs="?", const=2, help="perform UMAP on phenotypes and use resulting embeddings as phenotypes for GEMMA analysis \noptional: set amount of embeddings to be calculated (default: %(const)s) \nrecommended amount of embeddings: 1 - 5")
    parser.add_argument("-KC", "--kcpca", metavar="<float>", nargs='?', const=0.5, type=float, help="Kinship calculation via principal component analysis instead of GEMMA's internal method \noptional: r-squared threshold for LD pruning (default: %(const)s)")    
    parser.add_argument("-sv", "--sigval", metavar="<int>", type=int, default=6, help="set value where to draw significant line in manhattan plot \n<int> represents -log10(1e-<int>) (default: %(default)s) \nwhen using '-bslmm', value is adjusted to fit in the range of 0 to 1 \nset <int> to '0' to disable line")
    parser.add_argument("-fs", "--fontsize", metavar="<int>", type=int, default=26, help="Set the fontsize of plots \nDefault: %(default)s")
    parser.add_argument("-o", "--output", metavar="<path>", type=str, default=os.getcwd(), help="change the output directory \ndefault: %(default)s\ndirectory will be created if non-existent")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-lm", type=int, choices=[1,2,3,4], nargs="?", const=1, help="Association Tests with a Linear Model \noptional: specify which frequentist test to use (default: %(const)s) \n 1: performs Wald test \n 2: performs likelihood ratio test \n 3: performs score test \n 4: performs all three tests")
    group.add_argument("-gk", type=int, choices=[1,2], nargs="?", const=1, help="Estimate Relatedness Matrix from genotypes \noptional: specify which relatedness matrix to estimate (default: %(const)s) \n 1: calculates the centered relatedness matrix \n 2: calculates the standardized relatedness matrix")
    group.add_argument("-eigen", action="store_true", help="Perform Eigen-Decomposition of the Relatedness Matrix ")
    group.add_argument("-lmm", type=int, choices=[1,2,3,4], nargs="?", const=1, help="Association Tests with Univariate Linear Mixed Models \noptional: specify which frequentist test to use (default: %(const)s) \n 1: performs Wald test \n 2: performs likelihood ratio test \n 3: performs score test \n 4: performs all three tests \nTo perform Association Tests with Multivariate Linear Mixed Models, set '-multi' option")
    group.add_argument("-bslmm", type=int, choices=[1,2,3], nargs="?", const=1, help="Fit a Bayesian Sparse Linear Mixed Model \noptional: specify which model to fit (default: %(const)s) \n 1: fits a standard linear BSLMM \n 2: fits a ridge regression/GBLUP \n 3: fits a probit BSLMM")

    parser.add_argument("-ap", "--allphenotypes", action="store_true", help="all phenotypes will be used \nany phenotype selection with '-p' option will be overwritten")
    parser.add_argument("-ac", "--allcovariates", action="store_true", help="all covariates will be used \nany covariate selection with '-c' option will be overwritten")
    parser.add_argument("-m", "--multi", action="store_true", help="performs multivariate linear mixed model analysis with specified phenotypes \nonly active in combination with '-lmm' option")
    parser.add_argument("-nl", "--nolabel", action="store_true", help="remove the SNP labels in the manhattan plot \nreduces runtime if analysis results in many significant SNPs")
    parser.add_argument("-s", "--seed", action="store_true", help="perform UMAP with random seed \nreduces reproducibility")
    parser.add_argument("-r", "--retain", action="store_true", help="keep all temporary intermediate files \ne.g. subsetted and filtered VCF and .csv files")
    
    return parser.parse_args(argv)

class Parser:
    """Description:
    contains functions regarding argument parser.
    Parses arguments from command-line input.
    Class methods return variables from input arguments"""
    
    def __init__(self, args):
        self.args = getArgs(args)

    def set_geno(self):
        return self.args.vcf

    def set_pheno(self):
        return self.args.pfile

    def set_covar(self):
        return self.args.cfile

    def get_phenotypes(self):
        return self.args.pheno

    def get_covariates(self):
        return self.args.covar
    
    def set_gene_file(self):
        return self.args.genefile

    def set_gene_thresh(self):
        return self.args.genethresh

    def set_q(self):
        return self.args.minaf

    def set_n_top(self):
        return self.args.topsnp

    def set_memory(self):
        return self.args.memory

    def set_threads(self):
        return self.args.threads

    def set_A(self):
        return self.args.allphenotypes

    def set_B(self):
        return self.args.allcovariates

    def set_P(self):
        return self.args.PCA

    def set_U(self):
        return self.args.UMAP

    def set_sigval(self):
        return self.args.sigval

    def set_nolabel(self):
        return self.args.nolabel

    def set_fontsize(self):
        return self.args.fontsize

    def set_seed(self):
        return self.args.seed

    def set_keep(self):
        return self.args.retain

    def set_lm(self):
        return self.args.lm

    def set_gk(self):
        return self.args.gk

    def set_eigen(self):
        return self.args.eigen

    def set_lmm(self):
        return self .args.lmm

    def set_bslmm(self):
        return self.args.bslmm

    def set_filename(self):
        return self.args.relmatrix

    def set_pca(self):
        return self.args.kcpca

    def set_multi(self):
        return self.args.multi

    def set_out_dir(self):
        return self.args.output