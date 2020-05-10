configfile: "config.json"
BCFTOOLS_PLUGINS=config["BCFTOOLS_PLUGINS"]
plink2=config["plink2"]
eagle=config["eagle"]
Minimac3=config["Minimac3"]

rule all:
    input:
        expand("{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz_glm.PHENO1.glm.logistic.hybrid_manh.jpeg",chr=config["chromosome"],bfile=config["bfile"])

rule plink_to_vcf:
    input:
        bed="{bfile}.bed",
    	bim="{bfile}.bim",
    	fam="{bfile}.fam"
    params:
        prefix="{bfile}_chr{chr}"
    output:
        "{bfile}_chr{chr}.vcf"
    shell:
        '''
        plink2 --bed {input.bed} --bim {input.bim} --fam {input.fam} --chr {wildcards.chr} \
        --recode vcf --out {params.prefix}
        '''

rule wget_fasta:
    output:
        fasta="app/chr{chr}.fa.gz",
    params:
        prefix="app/chr{chr}.fa"
    shell:
        '''
        echo "Downloading fasta file, it may take a few minutes"
        wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr{wildcards.chr}.fa.gz -O {output.fasta}
        echo "Processing fasta file, it may take a few minutes"
        gunzip -f {output.fasta}
        sed -i "s/chr{wildcards.chr}/{wildcards.chr}/" {params.prefix}
        bgzip -f {params.prefix}
        samtools faidx {output.fasta}
        '''

rule wget_genetic_map_hg19:
    output:
        genetic_map_hg19="app/genetic_map_hg19.txt.gz",
    shell:
        '''
        echo "Downloading genetic_map_hg19 file, it may take a few minutes"
        wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19.txt.gz -O {output}
        '''

rule vcf_to_ref_vcf:
    input:
        fasta="app/chr{chr}.fa.gz",
        vcf="{bfile}_chr{chr}.vcf"
    output:
        "{bfile}_chr{chr}_ref.vcf.gz"
    priority: 50
    shell:
        '''
        export BCFTOOLS_PLUGINS={BCFTOOLS_PLUGINS}; \
        bcftools +fixref {input.vcf} \
        -- -f {input.fasta} -m flip -d | \
        bcftools sort -Ou | bcftools norm --rm-dup all -Oz -o {output}
        bcftools index {output}
        '''

rule wget_vcfRef:
    output:
        vcfRef="app/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        vcfRef_tbi="app/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
    shell:
        '''
        echo "Downloading reference panel vcf file, it may take a few minutes"
        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O {output.vcfRef}
        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi -O {output.vcfRef_tbi}
        '''

rule phasing:
    input:
        genetic_map_hg19="app/genetic_map_hg19.txt.gz",
        vcfTarget="{bfile}_chr{chr}_ref.vcf.gz",
        vcfRef="app/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        vcfRef_tbi="app/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
    output:
        "{bfile}_chr{chr}_ref_phased.vcf.gz"
    params:
        prefix="{bfile}_chr{chr}_ref_phased"
    shell:
        '''
        eagle \
        --geneticMapFile {input.genetic_map_hg19} \
        --vcfTarget {input.vcfTarget} \
        --vcfRef app/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --chrom {wildcards.chr} \
        --numThreads 1 \
        --outPrefix {params.prefix}
        '''

rule imputation:
    input:
        haps="{bfile}_chr{chr}_ref_phased.vcf.gz",
    output:
        "{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz"
    params:
        prefix="{bfile}_chr{chr}_ref_phased_imputed"
    shell:
        '''
        Minimac3 \
        --refHaps app/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --haps {input.haps} \
        --noPhoneHome --lowMemory \
        --prefix {params.prefix}
        '''

rule glm:
    input:
        "{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz"
    output:
        "{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz_glm.PHENO1.glm.logistic.hybrid"
    params:
        prefix="{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz_glm"
    shell:
        '''
        plink2 --vcf {input} dosage=DS \
        --id-delim _ \
        --glm sex hide-covar cols=chrom,pos,ref,alt,ax,test,nobs,orbeta,se,tz,p allow-no-covars \
        --update-sex {wildcards.bfile}.fam col-num=5 \
        --pheno <(awk '{{print $1,$2,$6}}' {wildcards.bfile}.fam) \
        --make-pfile \
        --out {params.prefix}
        '''

rule qqman:
    input:
        "{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz_glm.PHENO1.glm.logistic.hybrid"
    output:
        "{bfile}_chr{chr}_ref_phased_imputed.dose.vcf.gz_glm.PHENO1.glm.logistic.hybrid_manh.jpeg"
    script:
        "scripts/qqman.R"

### BK rules
#rule vcfRef_to_bcfRef:
#    input:
#        vcfRef="{vcfRef}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
#    output:
#        bcfRef="{vcfRef}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf",
#        bcfRef_csi="{vcfRef}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi"
#    shell: "bcftools view {input.vcfRef} -Ob -o {output.bcfRef}; bcftools index {output.bcfRef}"
