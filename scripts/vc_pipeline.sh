#!/usr/bin/env bash
# =============================================================================
# Per-sample short-read variant calling for pigeon pea (Cajanus cajan) WGS data
#
#   FastQC -> Trimmomatic -> BWA-MEM (+ read groups) -> sort -> MarkDuplicates
#          -> GATK HaplotypeCaller (GVCF mode)
#
# Usage:
#   bash scripts/vc_pipeline.sh config/config.sh              # all samples in the sheet
#   bash scripts/vc_pipeline.sh config/config.sh SRR12717516  # one sample only
#
# Samples whose GVCF already exists are skipped, so an interrupted run can simply be
# restarted. Then run scripts/joint_genotyping.sh to genotype and filter all samples together.
# =============================================================================
set -euo pipefail

CONFIG="${1:?Usage: bash scripts/vc_pipeline.sh <config.sh> [sample_id]}"
ONLY_SAMPLE="${2:-}"
# shellcheck source=/dev/null
source "$CONFIG"

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------- checks
for tool in $FASTQC $BWA $SAMTOOLS $GATK; do
    [[ "$tool" == java || "$tool" == -* || "$tool" == *.jar ]] && continue
    command -v "$tool" >/dev/null 2>&1 || fail "'$tool' not found. Activate the conda environment (environment.yml)."
done
[[ -f "$REF" ]]      || fail "Reference not found: $REF"
[[ -f "$SAMPLES" ]]  || fail "Sample sheet not found: $SAMPLES"
[[ -f "$ADAPTERS" ]] || fail "Adapter file not found: $ADAPTERS (see README, 'Adapters')"

mkdir -p "$OUTDIR"/{qc_reports,trimmed_fastq,aligned_bam,gvcf,logs}

# ---------------------------------------------------------------- step 0: reference indexes (once)
REF_BASE="${REF%.*}"
if [[ ! -f "${REF}.bwt" ]]; then
    log "[Step 0] Indexing reference for BWA (one-time)"
    $BWA index "$REF" > "$OUTDIR/logs/bwa_index.log" 2>&1
fi
[[ -f "${REF}.fai" ]] || $SAMTOOLS faidx "$REF"
if [[ ! -f "${REF_BASE}.dict" ]]; then
    log "[Step 0] Creating sequence dictionary for GATK (one-time)"
    $GATK --java-options "-Xmx${JAVA_MEM}" CreateSequenceDictionary -R "$REF" -O "${REF_BASE}.dict" \
        > "$OUTDIR/logs/create_dict.log" 2>&1
fi

# ---------------------------------------------------------------- per-sample function
process_sample() {
    local id="$1" r1="$2" r2="$3"
    local L="$OUTDIR/logs/$id"
    local gvcf="$OUTDIR/gvcf/${id}.g.vcf.gz"

    if [[ -s "$gvcf" && -s "${gvcf}.tbi" ]]; then
        log "[$id] GVCF already exists, skipping"
        return 0
    fi

    # Step 1: inputs present?
    [[ -f "$r1" && -f "$r2" ]] || fail "[$id] FASTQ files not found: $r1 / $r2"

    # Step 2: raw read quality
    log "[$id] Step 2/6  FastQC"
    $FASTQC -t "$THREADS" -o "$OUTDIR/qc_reports" "$r1" "$r2" > "${L}.fastqc.log" 2>&1

    # Step 3: adapter and quality trimming
    log "[$id] Step 3/6  Trimmomatic"
    local t="$OUTDIR/trimmed_fastq/$id"
    $TRIMMOMATIC PE -threads "$THREADS" "$r1" "$r2" \
        "${t}_1_paired.fastq.gz" "${t}_1_unpaired.fastq.gz" \
        "${t}_2_paired.fastq.gz" "${t}_2_unpaired.fastq.gz" \
        "ILLUMINACLIP:${ADAPTERS}:2:30:10" $TRIM_STEPS > "${L}.trimmomatic.log" 2>&1

    # Step 4: align, add read groups, sort (streamed: no SAM file is written)
    log "[$id] Step 4/6  BWA-MEM alignment + sort"
    local b="$OUTDIR/aligned_bam/$id"
    $BWA mem -t "$THREADS" -R "@RG\tID:${id}\tSM:${id}\tLB:${id}_lib1\tPL:ILLUMINA\tPU:${id}_unit1" \
        "$REF" "${t}_1_paired.fastq.gz" "${t}_2_paired.fastq.gz" 2> "${L}.bwa.log" \
      | $SAMTOOLS sort -@ "$THREADS" -o "${b}_sorted.bam" - 2> "${L}.sort.log"
    $SAMTOOLS index "${b}_sorted.bam"

    # Step 5: mark PCR/optical duplicates
    log "[$id] Step 5/6  MarkDuplicates"
    $GATK --java-options "-Xmx${JAVA_MEM}" MarkDuplicates \
        -I "${b}_sorted.bam" -O "${b}_dedup.bam" -M "${b}_dup_metrics.txt" \
        --CREATE_INDEX true > "${L}.markdup.log" 2>&1
    $SAMTOOLS flagstat -@ "$THREADS" "${b}_dedup.bam" > "${b}_flagstat.txt"

    # Step 6: per-sample variant calling in GVCF mode
    log "[$id] Step 6/6  GATK HaplotypeCaller (GVCF)"
    $GATK --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
        -R "$REF" -I "${b}_dedup.bam" -O "$gvcf" \
        -ERC GVCF --sample-ploidy "$PLOIDY" \
        --native-pair-hmm-threads "$THREADS" > "${L}.haplotypecaller.log" 2>&1

    if [[ "$KEEP_INTERMEDIATE" != true ]]; then
        rm -f "${t}"_*.fastq.gz "${b}_sorted.bam" "${b}_sorted.bam.bai"
    fi
    log "[$id] Done: $gvcf"
}

# ---------------------------------------------------------------- loop over the sample sheet
n=0
while IFS=$'\t' read -r id r1 r2 _; do
    [[ -z "$id" || "$id" == sample_id || "$id" == \#* ]] && continue
    [[ -n "$ONLY_SAMPLE" && "$id" != "$ONLY_SAMPLE" ]] && continue
    process_sample "$id" "$r1" "$r2"
    n=$((n + 1))
done < "$SAMPLES"

[[ $n -gt 0 ]] || fail "No samples processed (check $SAMPLES${ONLY_SAMPLE:+ for '$ONLY_SAMPLE'})"
log "Pipeline finished for $n sample(s). Next: bash scripts/joint_genotyping.sh $CONFIG"
