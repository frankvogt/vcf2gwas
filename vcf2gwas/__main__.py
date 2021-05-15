#!/usr/bin/env python

import sys
import os
import subprocess

from vcf2gwas.parsing import *

#from vcf2gwas.install import main as installer

argvals = None

def main(argvals=argvals):
    print("\nvcf2gwas v0.5 \n")
    print("Initialising..\n")
    #try:
    #    import umap
    #except:
    #    subprocess.run(["pip", "install", "umap-learn"])
    #if "umap-learn" not in str(subprocess.run(["conda", "list"], capture_output=True)):
    #    installer()
    string = str(subprocess.run(["conda", "list"], capture_output=True))
    modules = ["gemma", "plink", "bcftools"]
    if any(x not in string for x in modules):
        subprocess.run(["conda", "install", "-c", "bioconda", "bcftools==1.12*"])
        subprocess.run(["conda", "install", "-c", "bioconda", "plink==1.90*"])
        subprocess.run(["conda", "install", "-c", "bioconda", "gemma==0.98.3"])
    Parser(argvals)
    args = sys.argv[1:]
    #args.insert(0, 'conda')
    #args.insert(1, 'run')
    #args.insert(2, '--no-capture-output')
    #args.insert(3, '-n')
    #args.insert(4, 'vcf2gwas')
    #args.insert(5, 'python')
    #args.insert(6, os.path.join(os.path.dirname(__file__), 'starter.py'))
    args.insert(0, 'python')
    args.insert(1, os.path.join(os.path.dirname(__file__), 'starter.py'))
    subprocess.run(args)

if __name__ == '__main__':
    sys.exit(main())
