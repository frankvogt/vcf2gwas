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
import time

#from vcf2gwas.parsing import *
from parsing import *
from vcf2gwas.install import main as installer

argvals = None

def main(timestamp, argvals=argvals):
    
    version = set_version_number()
    print(f"\nvcf2gwas v{version} \n")
    print("Initialising..\n")
    P = Parser(argvals)
    args = sys.argv[1:]
    args.insert(0, 'python3.9')
    args.insert(1, os.path.join(os.path.dirname(__file__), 'starter.py'))
    
    try:
        args = delete_string(args, ['--timestamp'])
    except:
        pass
    args.insert(2, timestamp)
    args.insert(2, "--timestamp")
    
    geno = P.set_geno()
    if geno == "test":
        source = os.path.join(os.path.dirname(__file__), 'starter.py')
        installer()
        vcf = os.path.join("input", "example.vcf.gz")
        pheno = os.path.join("input", "example.csv")
        args = f'python3.9 {source} --timestamp {timestamp} -v {vcf} -pf {pheno} -p 1 -lm'.split()
    
    lm = P.set_lm()
    lmm = P.set_lmm()
    covar = P.set_covar()
    if lm == None and lmm == None:
        if covar != None:
            msg = "A covariate file can only be added when using the linear model ('-lm') or the linear mixed model ('-lmm')"
            raise SyntaxError(msg)

    process = subprocess.run(args)
    return process.returncode

def run_main():
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    try:
        sys.exit(main(timestamp))
    except KeyboardInterrupt as e:
        print("\nvcf2gwas interrupted")
        print("Cleaning up temporary files\n")
    finally:
        shutil.rmtree(f'_vcf2gwas_temp_{timestamp}', ignore_errors=True)

if __name__ == '__main__':
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    try:
        sys.exit(main(timestamp))
    except KeyboardInterrupt as e:
        print("\nvcf2gwas interrupted")
        print("Cleaning up temporary files\n")
    finally:
        shutil.rmtree(f'_vcf2gwas_temp_{timestamp}', ignore_errors=True)
