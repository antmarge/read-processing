#!/usr/bin/env python
import csv
import glob

# Where to output run results
RUN_DIR = "./"

# Where are the raw fastq files
RAW_DIR="./"

# A tab delimited file that lists 
# sampleID and sample filepath, one per line
sample_file = RAW_DIR + "/sample_list.txt"

# Dictionary of samples
from collections import defaultdict
sample_file_dict = defaultdict(list)

with open(sample_file) as f:
    for line in f:
       (key, val) = line.split()
       sample_file_dict[key].append(val)

SAMPLES = sample_file_dict.keys()

# Whether or not to trim read ends
# Not necessary (TRIM = 0) for modern DNA
# Trim = 2 for ancient DNA (partial UDG)
TRIM = ['0']

# LASER RUN ID
SEQ_RUN = "output-all"

# Directories

# ALL FASTQ FILES
FASTQ_DIR = RUN_DIR + "/1_fastq"
FASTQ_RAW_DIR = FASTQ_DIR + "/1a_raw"
FASTQ_ADAPTREM_DIR = FASTQ_DIR + "/1b_adaptrem"
FASTQ_TRIM_DIR = FASTQ_DIR + "/1c_trim"
TRIM_FASTQC_DIR = FASTQ_DIR + "/1d_fastqc"


# Alignment directories (BAMS)
BAM_DIR = RUN_DIR + "/2_bam"
BAM_ALN_DIR = BAM_DIR + "/2a_aligned"
BAM_MAPPED_DIR = BAM_DIR + "/2b_mapped"
BAM_RMDUP_DIR = BAM_DIR + "/2d_rmdup"
BAM_QUAL_DIR = BAM_DIR + "/2c_quality"

# Summary directories
COVERAGE_DIR = RUN_DIR + "/3_coverage"
SUMMARY_DIR = RUN_DIR + "/4_readSummary"

# Analysis directories
PSEUDO_DIR = RUN_DIR + "/5_pseudohaploid"
PCA_DIR = RUN_DIR + "/6_pca"
ADM_DIR = RUN_DIR + "/7_admixture"

# Human reference genome
REF = "/share/PI/pritch/Margaret/reference/hg19rCRS/hg19"

# Parameters
READ_MIN_LENGTH = "30"  # minimum read length after trimming 
BASE_MIN_QUAL = "30"  # base quality filter
OVERLAP = "1" # Minimum overlap between read and adapter allowed during Cutadapt
READ_START_TRIM = "2" #Number of bases to trim at beginning of all reads (minimizes C->T damaged bases)
READ_END_TRIM = "2" #Number of bases to trim at ends of all reads
SEED_DISABLE = "350"  # use a value bigger than read length to disable it
MAP_QUAL = "30"  # MAPQ filter
DEPTH = "17"  # Minimum number of reads covering a position to call variant
ADAPT_SEQ = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG" # 3' primer sequence for SE read

#### LASER read-based PCA
REF_VCF="/scratch/PI/pritch/AncientRome/analyses/modern_panels/preimputation/merged_matulloGlobe1kg_commonSNPs_global.vcf.gz"
DIR="/scratch/PI/pritch/AncientRome/analyses/laser/run5-global-20180530"
PILEUP2SEQ="/share/PI/pritch/Margaret/bin/LASER-2.04/pileup2seq/pileup2seq.py"
REF_FASTA="/share/PI/pritch/Margaret/reference/hg19g1k/human_g1k_v37.fa"
PANEL_NAME="matGlobe1kgGlobal"

# MapDamage parameters 
YMAX=0.15
DOWNSAMPLE=800000

#SOFTWARE
CUTADAPT="/share/PI/pritch/Margaret/bin/miniconda3/bin/cutadapt"
MAPDAMAGE="/home/users/antmarge/.local/bin/mapDamage"
MAPDAMAGE="/home/groups/pritch/Margaret/bin/mapDamage/mapDamage"
PICARD="/share/PI/pritch/Margaret/bin/picard-2.9.0/picard.jar"
DEPTH='/share/PI/pritch/Margaret/bin/depth-cover.jar'
SAMTOOLS = "/share/PI/pritch/Margaret/bin/samtools"

