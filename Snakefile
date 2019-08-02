# Snakefile for basic read processing from 
# raw fastq reads to alignment
# includes generation of read-based PCA (LASER) input files
# 
# This pipeline is appropriate for small files (< 2GB)
# There is a parallelized version where raw fastqs are split
# into 10 files each and processed in parallel
#
# Margaret Antonio
# Last updated: 2019.08.02

# Snakemake variables file specific to the run
include: "snakemake_variables.py"


import os,sys,re



rule all:
    input:
        expand(FASTQ_ADAPTREM_DIR + "/{sample}_adaptRem_trim{trim}.fastq.gz", sample = SAMPLES, trim = '0'),
        expand(BAM_MAPPED_DIR + "/{sample}_trim{trim}_mapped.bam", sample = SAMPLES, trim = '0'),
        expand(COVERAGE_DIR + "/{sample}/{sample}.summary.csv", sample = SAMPLES),
        expand(PCA_DIR + "/pileup/{sample}_trim0_mapped_q" + MAP_QUAL + "_rmdup.seq", sample = SAMPLES),
        expand(PCA_DIR + "/seq/" + SEQ_RUN + "_trim0_mapped_q" + MAP_QUAL + "_rmdup.seq")
	# Output files ehre
        
'''
rule fastqc_raw:
    input: lambda wildcards: sample_file_dict[{wildcards.sample}] 
    output: fastq_report = "{RUN_DIR}/" + FASTQC_FASTQ_RAW_DIR + "/{sample}_fastqc.html"
    benchmark: "{RUN_DIR}/benchmarks/fastqc_raw/{sample}.txt"
    shell:"""
        mkdir -p {RUN_DIR}
        mkdir -p {RUN_DIR}/reports
        mkdir -p {RUN_DIR}/reports/logs
        mkdir -p {RUN_DIR}/reports/reports
        mkdir -p {RUN_DIR}/{FASTQ_DIR}
        module load devel
        module load java/1.8.0_131
        mkdir -p {RUN_DIR}/{FASTQC_FASTQ_RAW_DIR}
        fastqc {input.fastq} --outdir={RUN_DIR}/{FASTQC_FASTQ_RAW_DIR}/ 
        """
'''


def getFastq(sample_id):
	return(sample_file_dict[sample_id])


# Local rules will not be submitted to the cluster
localrules: concat_fastqs


# Concatenate fastq files of same sample ID
rule concat_fastqs:
    input: lambda wildcards: sample_file_dict[wildcards.sample]
    output: fastq = FASTQ_RAW_DIR + "/{sample}.fastq.gz"
    benchmark: RUN_DIR + "/benchmarks/adaptRem/{sample}.txt"
    shell:"""
        cat {input} > {output}
    """

# Remove illumina universal adapters and trim reads if necessary
# Filter for min length
rule cutadapt:
    input: fastq = FASTQ_RAW_DIR + "/{sample}.fastq.gz"
    output: fastq = FASTQ_ADAPTREM_DIR + "/{sample}_adaptRem_trim{trim}.fastq.gz"
    benchmark: RUN_DIR + "/benchmarks/adaptRem/{sample}.txt"
    shell:"""
        {CUTADAPT} -u {wildcards.trim} -u -{wildcards.trim} -a {ADAPT_SEQ} \
            -m {READ_MIN_LENGTH} -O {OVERLAP} {input} | gzip -c > {output.fastq}
        """

# Run FASTQC on the fastq file that comes out of the cutadapt rule
# Good if uncertain about adapters are correct
rule fastqc_cutadapt:
    input: fastq = FASTQ_RAW_DIR + "/{sample}.fastq.gz"
    output: fastq_report = TRIM_FASTQC_DIR + "/trim{trim}/{sample}.adaptRem.trim{trim}.fastqc.html"
    benchmark: RUN_DIR + "benchmarks/fastqc/adaptRem_{sample}.txt"
    shell:"""
        module load devel
        module load java/1.8.0_131
        mkdir -p {FASTQC_FASTQ_ARTRIM_DIR}
        fastqc {input.fastq} --outdir={RUN_DIR}/{FASTQC_FASTQ_ARTRIM_DIR}/trim{wildcards.trim}
        echo $(date) ': Ran FASTQC on TRIMmed reads' >> {params.log}
        """

# Align reads with bwa aln
# BWA ALN used here since better for short reads (<50 bp)
# Use BWA MEM if reads are not short
# Then format for SAM using bwa samse
# Generate alignments in the SAM format given single-end reads
# Repetitive hits will be randomly chosen
# Only output mapped reads using -F4 flag

rule bwa_aln:
    input: FASTQ_ADAPTREM_DIR + "/{sample}_adaptRem_trim{trim}.fastq.gz"
    output: bam = BAM_MAPPED_DIR + "/{sample}_trim{trim}_mapped.bam"
    params:
        threads=5
    benchmark:  RUN_DIR + "/benchmarks/{sample}_trim{trim}_mapped.txt"
    shell:"""
        bwa aln -l {SEED_DISABLE} -t {params.threads} {REF} {input} | \
	bwa samse {REF} - {input} | {SAMTOOLS} view -Sbh - -F4 -@{params.threads} -b > {output.bam}
        """

# Sort reads, add read groups, filter for mapping quality
rule sort_rg_bam:
    input: BAM_MAPPED_DIR + "/{sample}_trim{trim}_mapped.bam"
    output: bam = BAM_QUAL_DIR + "/{sample}_trim{trim}_mapped_q" + MAP_QUAL + ".bam"
    benchmark: RUN_DIR + "/benchmarks/sortbam_{sample}.txt"
    shell:"""
        module load devel
        module load java/1.8.0_131
        module load python/3.6.1
        java -jar {PICARD} SortSam SORT_ORDER=coordinate COMPRESSION_LEVEL=0 \
                I={input} O=/dev/stdout VALIDATION_STRINGENCY=SILENT QUIET=TRUE | \
        java -jar {PICARD} AddOrReplaceReadGroups\
              I=/dev/stdin O=/dev/stdout VALIDATION_STRINGENCY=SILENT QUIET=TRUE \
              RGID={wildcards.sample} RGLB=liball COMPRESSION_LEVEL=0 RGPL=illumina RGPU=unit1 RGSM={wildcards.sample} | \
        {SAMTOOLS} view -h -q {MAP_QUAL} --output-fmt BAM - -o {output.bam} 
        {SAMTOOLS} index {output.bam}
        """

# Remove duplicates
rule rmdup:
    input: bam = BAM_QUAL_DIR + "/{sample}_trim{trim}_mapped_q" + MAP_QUAL + ".bam"
    output: bam = BAM_RMDUP_DIR + "/{sample}_trim{trim}_mapped_q" + MAP_QUAL + "_rmdup.bam"
    shell:"""
        {SAMTOOLS} rmdup --output-fmt BAM -s {input.bam} {output.bam}
        {SAMTOOLS} index {output.bam}
        """
# Get coverage
rule bam_coverage:
    input: bam = BAM_RMDUP_DIR + "/{sample}_trim0_mapped_q" + MAP_QUAL + "_rmdup.bam"
    output: 
        breakdown = COVERAGE_DIR + "/{sample}/{sample}.breakdown.csv",
        coverage = COVERAGE_DIR + "/{sample}/{sample}.coverage.csv",
        summary = COVERAGE_DIR + "/{sample}/{sample}.summary.csv"
    params:
        dir = COVERAGE_DIR + "/{sample}",
    shell:"""
        module load devel
        module load java/1.8.0_131
        mkdir -p {params.dir}
        java -jar {DEPTH} -bam {input.bam} -n {wildcards.sample} -d {params.dir}
        """


# Pileup for read-based PCA
rule mpileup:
    input: bam = BAM_RMDUP_DIR + "/{sample}_trim0_mapped_q" + MAP_QUAL + "_rmdup.bam",
           bed = "/share/PI/pritch/Margaret/bin/LASER-2.04/HGDP/HGDP_938_withchr.bed"
    output: PCA_DIR + "/pileup/{sample}_trim0_mapped_q" + MAP_QUAL + "_rmdup.pileup"
    shell: """
        {SAMTOOLS} mpileup -q 30 -Q 20 -f {REF}.fa -l {input.bed} {input.bam}  > {output}
    """

# MERGE Sample individuals mpileups to get seq and site files
# The pileup2seq.py script from LASER requires python 2


# Creates the seq file for input into LASER
rule pileup2seq:
    input: pileups = expand(PCA_DIR + "/pileup/{sample}_trim0_mapped_q" + MAP_QUAL + "_rmdup.pileup",sample = SAMPLES),
           site = "/share/PI/pritch/Margaret/bin/LASER-2.04/HGDP/HGDP_938_withchr.site"
    output: PCA_DIR + "/seq/" + SEQ_RUN + "_trim0_mapped_q" + MAP_QUAL + "_rmdup.seq"
    params: outroot = lambda wildcards, output: output[0][:-4]
    shell: """
	source activate python27
	python {PILEUP2SEQ} -m {input.site} -f {REF}.fa -o {params.outroot} {input.pileups}
	source deactivate
    """

# Run laser with the example HGDP reference panel included in laser
