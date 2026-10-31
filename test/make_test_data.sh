#!/usr/bin/env bash
# Build a small simulated dataset (3 samples, ~40x coverage, 250 kb genome) in test/data/.
# Needs: python3, wgsim (bioconda package "wgsim", or bundled with older samtools).
set -euo pipefail
cd "$(dirname "$0")"

python3 simulate_population.py

for s in geno_01 geno_02 geno_03; do
    # ~33,000 read pairs of 2 x 150 bp from the two haplotypes of each sample,
    # 0.5% sequencing error, no extra random mutations (-r 0 -R 0)
    seed=$(( 10#${s#geno_} * 7 ))
    wgsim -N 33000 -1 150 -2 150 -e 0.005 -r 0 -R 0 -X 0 -S "$seed" \
          "data/${s}_hap.fa" "data/${s}_1.fastq" "data/${s}_2.fastq" > /dev/null
    gzip -f "data/${s}_1.fastq" "data/${s}_2.fastq"
done

printf "sample_id\tfastq_1\tfastq_2\n" > data/samples_test.tsv
for s in geno_01 geno_02 geno_03; do
    printf "%s\ttest/data/%s_1.fastq.gz\ttest/data/%s_2.fastq.gz\n" "$s" "$s" "$s" >> data/samples_test.tsv
done
echo "Test data ready in test/data/"
