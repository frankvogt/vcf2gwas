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

import shutil
import sys
import os
import subprocess

from vcf2gwas.parsing import *
from vcf2gwas.install import main as installer

argvals = None

def main(argvals=argvals):
    
    print("\nvcf2gwas v0.7.4 \n")
    print("Initialising..\n")
    P = Parser(argvals)
    args = sys.argv[1:]
    args.insert(0, 'python3.9')
    args.insert(1, os.path.join(os.path.dirname(__file__), 'starter.py'))
    
    geno = P.set_geno()
    if geno == "test":
        source = os.path.join(os.path.dirname(__file__), 'starter.py')
        installer()
        vcf = os.path.join("input", "example.vcf.gz")
        pheno = os.path.join("input", "example.csv")
        args = f'python3.9 {source} -v {vcf} -pf {pheno} -p 1 -lm'.split()
    
    lm = P.set_lm()
    lmm = P.set_lmm()
    covar = P.set_covar()
    if lm == None and lmm == None:
        if covar != None:
            sys.exit(print("Error: A covariate file can only be added when using the linear model ('-lm') or the linear mixed model ('-lmm')"))

    subprocess.run(args)

    shutil.rmtree("_vcf2gwas_temp", ignore_errors=True)

if __name__ == '__main__':
    sys.exit(main())
