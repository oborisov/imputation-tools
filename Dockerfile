FROM ubuntu:18.04
FROM continuumio/miniconda3

# installing basic command line tools and libs
RUN apt-get install -y wget tar libgomp1

# downloading and extracting data
RUN wget https://github.com/oborisov/imputation-tools/raw/master/Snakefile && \
wget https://github.com/oborisov/imputation-tools/raw/master/config.json && \
wget https://github.com/oborisov/imputation-tools/raw/master/plink2_eagle_Minimac3.tar.gz -P /app && \
wget https://github.com/oborisov/imputation-tools/raw/master/simulated_test_data.tar.gz -P /data && \
tar -xzvf /app/plink2_eagle_Minimac3.tar.gz -C /app && mv /app/plink2_eagle_Minimac3/* /app/ && \
tar -xzvf /data/simulated_test_data.tar.gz -C /data && mv /data/simulated_test_data/* /data

# installing conda software
RUN conda create -c conda-forge -c bioconda -n snakemake snakemake
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]
RUN conda install -c bioconda bcftools samtools tabix

# installing r and packages
RUN conda install -c r -c bioconda r-base r-data.table r-qqman

# setting binaries as executable
RUN chmod 755 /app/Minimac3 /app/eagle /app/plink2
ENV PATH "$PATH:/app"

# running imputation-tools
#RUN PATH=$PATH:/app && snakemake
