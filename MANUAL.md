# vcf2gwas manual

## About

This manual is meant to clarify the most important options and possibilities when running vcf2gwas.  
To obtain example files for running vcf2gwas, run:

```
vcf2gwas -v test
```

Besides testing the installation, this will copy the example VCF file `example.vcf.gz` and phenotype file `example.csv` to your current working directory.  
The VCF file `example.vcf.gz` contains pre-filtered SNP information of about 60 accessions of *A. thaliana* from the [1001 Genomes](https://1001genomes.org/accessions.html) project.  
The phenotype file `example.csv` is made up of the `avrRpm1` phenotype from the [AraPheno](https://arapheno.1001genomes.org/phenotype/17/) database for the same accessions.

## Contents

* [GEMMA models](#gemma-models)
    * [Running a linear model analysis](#running-a-linear-model-analysis)
    * [Running a linear mixed model analysis](#running-a-linear-mixed-model-analysis)
    * [Running a bayesian sparse linear mixed model analysis](#running-a-bayesian-sparse-linear-mixed-model-analysis)
* [File related options](#file-related-options)
    * [Selecting multiple phenotype files](#selecting-multiple-phenotype-files)
    * [Analyzing multiple phenotypes](#analyzing-multiple-phenotypes)
    * [Adding covariates](#adding-covariates)
    * [Comparing results to specific genes](#comparing-results-to-specific-genes)
    * [Adding a relatedness matrix](#adding-a-relatedness-matrix)
    * [Changing the output directory](#changing-the-output-directory)
* [Miscellaneous options](#miscellaneous-options)
    * [Limiting memory and core usage](#limiting-memory-and-core-usage)
    * [Using dimensionality reduction of phenotypes for analysis](#using-dimensionality-reduction-of-phenotypes-for-analysis)
    * [Using PCA of genotypes instead of standard relatedness matrix](#using-pca-of-genotypes-instead-of-standard-relatedness-matrix)
    * [Filter out SNPs](#filter-out-snps)
    * [Change manhattan plot](#change-manhattan-plot)
    * [Change font size in all plots](#change-font-size-in-all-plots)
    * [Keep temporary files](#keep-temporary-files)
* [Output](#output)
    * [Plots](#plots)
    * [Summaries](#summaries)

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

vcf2gwas will summarize the `n` best SNPs (specified with `-t/--topsnp`) of every analyzed phenotype and compare them to the genes in the file by calculating the distance between each SNP and gene upstream as well as downstream. These results can be filtered by saving only those SNPs with a distance to a gene lower than a specific threshold (set with `-gt/--genethresh`).

### Adding a relatedness matrix

Although vcf2gwas will by default calculate a relatedness matrix depending on the chosen model, one may want to add a different one instead to the analysis. This is possible by employing the `-k/--relmatrix` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -k [filename]
```

To use vcf2gwas to just calculate a relatedness matrix from the VCF file, run the `-gk` option:

```
vcf2gwas -v [filename] -gk 
```

To calculate the relatedness matrix and perform its eigen-decomposition in the same run, use the `-eigen` option:

```
vcf2gwas -v [filename] -eigen
```

Of course the `-eigen` option can also be used when supplying your own relatedness matrix with the `-k/--relmatrix` option.

### Changing the output directory

By default, vcf2gwas will save the output in the current working directory. To change to a unique output directory, use the `-o/--output` option to specify a path:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -o dir/example/
```

## Miscellaneous options

In this section other useful options of vcf2gwas will be elucidated.

### Limiting memory and core usage

By default vcf2gwas will use half of the available memory and all logical cores minus one. It can be important to limit usage of these resources especially when running the analysis on a machine shared with others.  
To set the memory (in MB) and core usage employ the `-M/--memory` and `-T/--threads` option, respectively:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -M 8000 -T 6
```

Now, vcf2gwas uses 8 GB of memory and 6 cores to carry out the analysis.  
It is recommended to not set the memory to less than 1 GB.

### Using dimensionality reduction of phenotypes for analysis

When analyzing many phenotypes it can be escpecially beneficial to reduce the phenotypic dimensions. This allows the user to analyze any underlying structure in their phenotypic data by using the output of the dimensionality reduction as phenotypes for GEMMA.  
vcf2gwas offers to often-used methods to reduce the dimensions: principal component analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)) and Uniform Manifold Approximation and Projection ([UMAP](https://arxiv.org/abs/1802.03426)).  
Both can be used either separately or simultaneously in the analysis.

#### PCA

To perform PCA on the phenotype data and use the principal components as phenotypes for the analysis, use the `-P/--PCA` option.  
By default, vcf2gwas will reduce the phenotype dimensionality to 2 PCs. To change this value to any value between 2 and 10, append the value to the option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -P 3
```

Now, vcf2gwas will reduce the phenotype dimensionality to 3 instead of 2.

#### UMAP

To perform UMAP reduction on the phenotype data and use the embeddings as phenotypes for the analysis, use the `-U/--UMAP` option.  
By default, vcf2gwas will reduce the phenotype dimensionality to 2 embeddings. To change this value to any value between 1 and 5, append the value to the option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -U 3
```

Now, vcf2gwas will reduce the phenotype dimensionality to 3 instead of 2.

### Using PCA of genotypes instead of standard relatedness matrix

vcf2gwas uses GEMMAs standard method of kinship calculation for the linear mixed model, which produces a relatedness matrix. Instead of using this standard method, the relatedness matrix can optionally be calculated via PCA by utilizing the `-KC/--kcpca` option.  
The SNP data from the VCF file will be pruned by linkage disequilibrium with a default r-squared threshold of 0.5. To change the threshold, append the value to the option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -KC 0.8
```

### Filter out SNPs

By default, vcf2gwas will filter out SNPs with a minimum allele frequency of 0.01. To change this threshold use the `-q/--minaf` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -q 0.02
```

### Change manhattan plot

The manhattan plot which will be produced from the GEMMA output, has by default a significant value line drawn at *-log10(1e-6)*. All SNPs above that line will be labeled.  
To change this threshold line, use the `-sv/--sigval` option.

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -sv 7
```

The line will now be drawn at *-log10(1e-7)*.  
To disable the line and not label any SNPs, change the value to 0.

To remove the SNP lables completely, utilize the `-nl/--nolabel` option.  

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -nl
```

This can be beneficial to reduce the overall runtime when the analysis results in many significant SNPs.  

### Change font size in all plots

Since the resulting plots may be used in various contexts, the font size of the plots produced by vcf2gwas can be changed by using the `-fs/--fontsize option.

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -fs 20
```

The font sizes of the plots will now be changed from the default value of 26 pt to 20 pt.

### Keep temporary files

During the analysis, various temporary files like subsetted and filtered VCF and `.csv` files are produced. By default they are removed once they are no longer needed but if one wants to retain these files, employ the `-r/--retain` option:

```
vcf2gwas -v [filename] -pf [filename] -p 1 -lmm -r
```

## Output

The following part shows some of the output plots and summaries produced by the analysis of the example files using the linear mixed model and standard options.

### Plots

Manhattan-plot labeling significant SNPs with a standard significant value threshold of *-log10(1e-6)*:
<img src="https://github.com/frankvogt/vcf2gwas/blob/main/files/lmm_manh.png" alt="Manhattan-plot" width="75%"/>

QQ-plot comparing the expected and observed probability distributions:
<img src="https://github.com/frankvogt/vcf2gwas/blob/main/files/lmm_qq.png" alt="QQ-plot" width="75%"/>

### Summaries

Amongst other things, vcf2gwas will sort the SNPs of every analyzed phenotype, save the specified amount of top SNPs for each phenotype, summarize these SNPs of all phenotypes to check if certain SNPs occur more than once and optionally compare these SNPs to the genes supplied by the gene file. Below is the [output](https://github.com/frankvogt/vcf2gwas/blob/main/files/compared_summarized_top_SNPs_complete_example.csv) shown of such a gene comparison when supplying a [gene file](https://github.com/frankvogt/vcf2gwas/blob/main/files/NLR.csv) containing information about NLR genes in *A. thaliana*:

|SNP_ID|chr|phenotypes|gene_ID(up)|gene_comment(up)|gene_name(up)|gene_distance(up)|SNP_pos|gene_distance(down)|gene_name(down)|gene_comment(down)|gene_ID(down)|
|---|---|---|---|---|---|---|---|---|---|---|---|
|3:2237364|3|avrRpm|AT3G07040.1|NB-ARC domain-containing disease resistance protein|RPM1|8340|2237364|||||
|3:2237394|3|avrRpm|AT3G07040.1|NB-ARC domain-containing disease resistance protein|RPM1|8370|2237394|||||
|3:2237446|3|avrRpm|AT3G07040.1|NB-ARC domain-containing disease resistance protein|RPM1|8422|2237446|||||
|3:2237452|3|avrRpm|AT3G07040.1|NB-ARC domain-containing disease resistance protein|RPM1|8428|2237452|||||
|3:2289171|3|avrRpm|AT3G07040.1|NB-ARC domain-containing disease resistance protein|RPM1|60147|2289171|||||
|2:975138|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||84618|975138|28330||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|2:975234|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||84714|975234|28234||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|2:975320|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||84800|975320|28148||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|2:975405|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||84885|975405|28063||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|2:975411|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||84891|975411|28057||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|2:975490|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||84970|975490|27978||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|2:975590|2|avrRpm|AT2G03030.1|Toll-Interleukin-Resistance (TIR) domain family protein||85070|975590|27878||Toll-Interleukin-Resistance (TIR) domain family protein|AT2G03300.1|
|1:11273854|1|avrRpm|||||11273854|14598|RAC1|Disease resistance protein (TIR-NBS-LRR class) family|AT1G31540.2|
|1:11273813|1|avrRpm|||||11273813|14639|RAC1|Disease resistance protein (TIR-NBS-LRR class) family|AT1G31540.2|
|3:2165688|3|avrRpm|||||2165688|60264|RPM1|NB-ARC domain-containing disease resistance protein|AT3G07040.1|

