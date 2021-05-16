# vcf2gwas manual

## About

This manual is meant to clarify the most important options and possibilities when running vcf2gwas.  
To obtain example files for running vcf2gwas, run:

```
vcf2gwas -v test
```

Besides testing the installation, this will copy the example `VCF` file `example.vcf.gz` and phenotype file `example.csv` to your current working directory.


## GEMMA models

This sections demonstrates the usage of vcf2gwas when running the different analysis models of GEMMA.  

To read more in detail about the different models GEMMA is able to run, please refer to their [manual](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.pdf).

### Running a linear model analysis

To run a linear model analysis with default settings, select a `VCF` file with the `-v/--vcf` and a phenotype file with the `-pf/--pfile` option, respectively. To select a specific phenotype from the phenotype file, utilize the `-p/--pheno` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lm
```

This will use the first phenotype column in the phenotype file for analysis.

### Running a linear mixed model analysis

GEMMA offers univariate and multivariate linear mixed model analysis. 

#### Univariate linear mixed model

To run a linear mixed model analysis with default settings, type:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm
```

This again uses the first phenotype column in the phenotype file to carry out the analysis.  
To analyse multiple phenotypes (for example the 1st, 2nd and 4th) independently, run:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -p 2 -p 4 -lmm
```

This will carry out the univariate linear mixed model analysis for the 1st, 2nd and 4th phenotype column in the phenotype file.

#### Multivariate linear mixed model

To analyse these phenotypes jointly, employ the `-m/--multi` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -p 2 -p 4 -lmm --multi
```

### Running a bayesian sparse linear mixed model analysis

Run a bayesian sparse linear mixed model analysis with any phenotype in the phenotype file:

```
vcf2gwas -v [filename] -pf [filename] -p [int] -bslmm
```

## File related options

vcf2gwas is capable of analyzing multiple phenotypes from one or multiple phenotypes at once. 
Depending on the available cores and memory of the machine running the analysis and the number of phenotypes to be analyzed, the phenotype file will be split up during the analysis to ensure maximum efficiency by analyzing phenotypes in parallel.

### Selecting multiple phenotype files

vcf2gwas is able to take in multiple phenotype files by employing the `-pf/--pfile` option for each file:

```
vcf2gwas -v [filename] -pf [filename1] -pf [filename2] -p [int] -lmm
```

If multiple phenotypes are specified, these phenotypes will be analyzed in every phenotype file.
This is advantageous if one wants to analyze the same phenotypes with a different subset of the individuals in the `VCF` file. 

### Analyzing multiple phenotypes

To analyze multiple phenotypes in one run you can either specify multiple phenotypes

```
vcf2gwas -v [filename] -pf [filename] -p [num1] -p [num2] -p [num3] -lmm
```

or select all phenotypes in the phenotype file at once utilizing the `-ap/--allphenotypes` option:

```
vcf2gwas -v [filename] -pf [filename] -ap -lmm
```

### Adding covariates

GEMMA supports adding covariates to the linear mixed model analysis. 
Selecting covariates from a covariate file follows the same scheme as selecting phenotypes, using `-c/--covar` and `-cf/--cfile`, respectively:

```
vcf2gwas -v [filename] -pf [filename] -p [num1] -cf [filename] -c 1 -lmm
```

Here, the 1st covariate column of the covariate file will be considered in the analysis of the selected phenotype.
Similarly to the phenotype options, multiple covariates can be selected as well as all at once using `-ac/--allcovariates`.
If covariates are selected when running a different GEMMA model, the covariates will simply be added to the analysis as phenotypes. 

### Comparing results to specific genes

Comparing the GWAS results to specific genes of interest can be a tedious task. To facilitate this process, vcf2gwas supports adding a 'gene-file' containing the position as well as additional information of genes to the analysis by using the `-gf/--genefile` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -gf [filename]
```

vcf2gwas will summarize the `n` best SNPs (specified with `-t/--topsnp`) of every analyzed phenotype and compare them to the genes in the file by calculating the distance between each SNP and gene. These results can be filtered by saving only those SNPs with a distance to a gene lower than a specific threshold (set with `-gt/--genethresh`).

### Add relatedness matrix

Although vcf2gwas will by default calculate a relatedness matrix depending on the chosen model, one may want to add a different one to the analysis. This is possible by employing the `-k/--relmatrix` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -k [filename]
```

### Change the output directory

By default, vcf2gwas will save the output in the current working directory. To change to a unique output directory, use the `-o/--output` option to specify a path:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -o dir/example/
```

## Miscellaneous options
