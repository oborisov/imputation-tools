# imputation-tools
A set of tools to impute high throughput genotyping data. Data in [plink binary format](https://www.cog-genomics.org/plink/1.9/formats#bed) is expected as input. *imputation-tools* will:
1. Split the data by chromosomes
2. Align the data according to reference strand using [bcftools fixref](https://samtools.github.io/bcftools/howtos/plugin.fixref.html)
3. Phase the data with [Eagle2](https://data.broadinstitute.org/alkesgroup/Eagle/) and impute the data using [Minimac3](https://genome.sph.umich.edu/wiki/Minimac3) using [1000 genomes reference panel](https://data.broadinstitute.org/alkesgroup/Eagle/#x1-300005.3)
4. Perform standard association test with [plink2 --glm](https://www.cog-genomics.org/plink/2.0/assoc#glm)
5. Visualize the GWAS results using [qqman](https://cran.r-project.org/web/packages/qqman/index.html)

## installation
*imputation-tools* was tested on Ubuntu 18.04.4 LTS. Using 2.50GHz CPU (1 thread) the whole analysis for chromosome 22 takes approximately 90 minutes. *imputation-tools* requires packages that can be installed using conda ([a short guide to install conda](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1-installing-miniconda-3))

clone repoitory, unpack files, install environment and required packages, activate environment:
```
git clone https://github.com/oborisov/imputation-tools.git  
cd imputation-tools
gunzip data/* app/*
conda env create --file imputation-tools.yml  
conda activate imputation-tools
```

## Running
```
/home/user/miniconda3/envs/imputation-tools/libexec/bcftools/
cd imputation-tools
snakemake --config chromosome=22 --config bfile=data/sim1_GSA --config BCFTOOLS_PLUGINS=$(which bcftools | sed 's/bin\/bcftools/libexec\/bcftools/')
```
The following keys set the options for the *imputation-tools*:
chromosome to be imputed: ```--config chromosome=22```  
path to [binary plink file](https://samtools.github.io/bcftools/howtos/plugin.fixref.html) prefix: ```--config bfile=data/sim1_GSA```  
path to [BCFTOOLS_PLUGINS](https://samtools.github.io/bcftools/howtos/plugins.html), should be determined automatically based on the bcftools installation via conda: ```--config BCFTOOLS_PLUGINS=$(which bcftools | sed 's/bin\/bcftools/libexec\/bcftools/')```  

### used software
snakemake: ```snakemake.readthedocs.io/```  
plink2: ```https://www.cog-genomics.org/plink/2.0/```  
bcftools: ```http://samtools.github.io/bcftools/bcftools.html```  
samtools: ```http://samtools.github.io/```  
tabix: ```http://www.htslib.org/doc/tabix.html```  
Eagle2: ```https://data.broadinstitute.org/alkesgroup/Eagle/```  
Minimac3: ```https://genome.sph.umich.edu/wiki/Minimac3```  


## Docker (under development)
pull latest version of container: ```docker pull olegborisov/imputation-tools:latest```  
run application: ```docker run olegborisov/imputation-tools conda run -n snakemake snakemake```
