#!/bin/bash --login
set -e

conda activate python-app

exec "vcf2gwas" "$@"
