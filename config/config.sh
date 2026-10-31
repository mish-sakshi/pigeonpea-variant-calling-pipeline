# shellcheck shell=bash disable=SC2034
# =============================================================================
# Pipeline configuration. Edit this file only; the scripts read everything from here.
# Paths are relative to the folder you run the pipeline from (the repo root).
# =============================================================================

# ---- inputs ----
REF="reference/pigeonpea_reference.fa"   # Cajanus cajan reference genome (FASTA)
SAMPLES="config/samples.tsv"             # sample sheet: sample_id <tab> fastq_1 <tab> fastq_2

# ---- outputs ----
OUTDIR="results"                         # everything the pipeline writes goes here
COHORT_NAME="pigeonpea_cohort"           # prefix for the joint-genotyped VCF

# ---- compute ----
THREADS=16
JAVA_MEM="8g"                            # memory for GATK / Picard (Java heap)

# ---- trimming (Trimmomatic) ----
ADAPTERS="adapters/TruSeq3-PE.fa"        # adapter file shipped with Trimmomatic
TRIM_STEPS="SLIDINGWINDOW:4:20 MINLEN:50"

# ---- variant calling ----
PLOIDY=2

# ---- tool commands (change if your installation names them differently) ----
FASTQC="${FASTQC:-fastqc}"
TRIMMOMATIC="${TRIMMOMATIC:-trimmomatic}"                # e.g. "java -jar /path/to/trimmomatic-0.39.jar"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
GATK="${GATK:-gatk}"

# ---- housekeeping ----
KEEP_INTERMEDIATE=false                  # true = keep trimmed FASTQs and unmarked BAMs
