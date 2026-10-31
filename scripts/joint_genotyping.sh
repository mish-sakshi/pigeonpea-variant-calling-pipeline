#!/usr/bin/env bash
# =============================================================================
# Joint genotyping and hard filtering of all per-sample GVCFs
#
#   CombineGVCFs -> GenotypeGVCFs -> split SNPs / INDELs -> GATK hard filters -> PASS-only VCF
#
# Usage:  bash scripts/joint_genotyping.sh config/config.sh
# Filters follow GATK's generic hard-filtering recommendations for organisms without a
# truth/training variant set (no VQSR possible for pigeon pea).
# =============================================================================
set -euo pipefail

CONFIG="${1:?Usage: bash scripts/joint_genotyping.sh <config.sh>}"
# shellcheck source=/dev/null
source "$CONFIG"

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }
G() { $GATK --java-options "-Xmx${JAVA_MEM}" "$@"; }

V="$OUTDIR/variants"
L="$OUTDIR/logs/joint"
mkdir -p "$V"

# Collect GVCFs in sample-sheet order
args=()
while IFS=$'\t' read -r id _; do
    [[ -z "$id" || "$id" == sample_id || "$id" == \#* ]] && continue
    g="$OUTDIR/gvcf/${id}.g.vcf.gz"
    [[ -s "$g" ]] || fail "Missing GVCF for $id ($g). Run scripts/vc_pipeline.sh first."
    args+=(--variant "$g")
done < "$SAMPLES"
log "Joint genotyping $(( ${#args[@]} / 2 )) sample(s)"

log "Step 1/4  CombineGVCFs"
G CombineGVCFs -R "$REF" "${args[@]}" -O "$V/${COHORT_NAME}.g.vcf.gz" > "${L}.combine.log" 2>&1

log "Step 2/4  GenotypeGVCFs"
G GenotypeGVCFs -R "$REF" -V "$V/${COHORT_NAME}.g.vcf.gz" -O "$V/${COHORT_NAME}.raw.vcf.gz" \
    > "${L}.genotype.log" 2>&1

log "Step 3/4  Hard filtering (SNPs and INDELs separately)"
G SelectVariants -V "$V/${COHORT_NAME}.raw.vcf.gz" --select-type-to-include SNP \
    -O "$V/${COHORT_NAME}.snps.vcf.gz" > "${L}.select_snps.log" 2>&1
G VariantFiltration -V "$V/${COHORT_NAME}.snps.vcf.gz" -O "$V/${COHORT_NAME}.snps.filtered.vcf.gz" \
    -filter "QD < 2.0"              --filter-name "QD2" \
    -filter "QUAL < 30.0"           --filter-name "QUAL30" \
    -filter "SOR > 3.0"             --filter-name "SOR3" \
    -filter "FS > 60.0"             --filter-name "FS60" \
    -filter "MQ < 40.0"             --filter-name "MQ40" \
    -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    > "${L}.filter_snps.log" 2>&1

G SelectVariants -V "$V/${COHORT_NAME}.raw.vcf.gz" --select-type-to-include INDEL \
    -O "$V/${COHORT_NAME}.indels.vcf.gz" > "${L}.select_indels.log" 2>&1
G VariantFiltration -V "$V/${COHORT_NAME}.indels.vcf.gz" -O "$V/${COHORT_NAME}.indels.filtered.vcf.gz" \
    -filter "QD < 2.0"               --filter-name "QD2" \
    -filter "QUAL < 30.0"            --filter-name "QUAL30" \
    -filter "FS > 200.0"             --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    > "${L}.filter_indels.log" 2>&1

log "Step 4/4  Merge and keep PASS variants"
G MergeVcfs -I "$V/${COHORT_NAME}.snps.filtered.vcf.gz" -I "$V/${COHORT_NAME}.indels.filtered.vcf.gz" \
    -O "$V/${COHORT_NAME}.filtered.vcf.gz" > "${L}.merge.log" 2>&1
G SelectVariants -V "$V/${COHORT_NAME}.filtered.vcf.gz" --exclude-filtered \
    -O "$V/${COHORT_NAME}.PASS.vcf.gz" > "${L}.select_pass.log" 2>&1

# Small summary table
summary="$V/${COHORT_NAME}.summary.tsv"
{
    printf "set\tSNPs\tINDELs\n"
    for f in raw PASS; do
        zcat "$V/${COHORT_NAME}.${f}.vcf.gz" | awk -v set="$f" '
            !/^#/ { snp = (length($4) == 1); n = split($5, alt, ",")
                    for (k = 1; k <= n; k++) if (length(alt[k]) != 1) snp = 0
                    if (snp) s++; else i++ }
            END   { printf "%s\t%d\t%d\n", set, s, i }'
    done
} > "$summary"
cat "$summary"
log "Done. Final call set: $V/${COHORT_NAME}.PASS.vcf.gz"
