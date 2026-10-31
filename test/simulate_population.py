#!/usr/bin/env python3
"""
Simulate a tiny diploid "population" for testing the pipeline.

Writes to test/data/:
  ref.fa                      random 2-contig reference (150 kb + 100 kb)
  <sample>_hap.fa             the sample's two haplotypes (used by wgsim to make reads)
  truth.tsv                   every simulated variant with each sample's genotype
  TruSeq3-PE.fa               Trimmomatic adapter file
No real pigeon pea data is used.
"""
import os
import random

random.seed(2025)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)

CONTIGS = {"chr_test1": 150_000, "chr_test2": 100_000}
SAMPLES = ["geno_01", "geno_02", "geno_03"]
N_SITES = {"chr_test1": 180, "chr_test2": 120}     # variant sites in the population
INDEL_FRACTION = 0.15

ref = {c: [random.choice("ACGT") for _ in range(n)] for c, n in CONTIGS.items()}

# Choose variant sites at least 30 bp apart; build REF/ALT alleles
sites = []
for c, n in N_SITES.items():
    positions = sorted(random.sample(range(200, CONTIGS[c] - 200, 30), n))
    for pos in positions:                                      # 0-based position
        r = ref[c][pos]
        if random.random() < INDEL_FRACTION:
            L = random.randint(1, 5)
            if random.random() < 0.5:                          # insertion
                alt = r + "".join(random.choice("ACGT") for _ in range(L))
                refallele = r
            else:                                              # deletion
                refallele = "".join(ref[c][pos:pos + L + 1])
                alt = r
        else:
            refallele, alt = r, random.choice([b for b in "ACGT" if b != r])
        # genotype per sample: 0 = hom-ref, 1 = het, 2 = hom-alt; every site varies in >=1 sample
        while True:
            gts = [random.choices([0, 1, 2], weights=[0.5, 0.3, 0.2])[0] for _ in SAMPLES]
            if any(gts):
                break
        sites.append((c, pos, refallele, alt, gts))


def haplotype(contig, sample_idx, hap):
    """Apply this sample's alleles to one haplotype (hap 0 or 1)."""
    seq, out, i = ref[contig], [], 0
    for c, pos, ra, alt, gts in sites:
        if c != contig:
            continue
        g = gts[sample_idx]
        carries = (g == 2) or (g == 1 and hap == 0)          # het alleles on haplotype 0
        if carries:
            out.extend(seq[i:pos])
            out.append(alt)
            i = pos + len(ra)
    out.extend(seq[i:])
    return "".join(out)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, s in records:
            fh.write(f">{name}\n")
            for k in range(0, len(s), 70):
                fh.write(s[k:k + 70] + "\n")


write_fasta(f"{OUT}/ref.fa", [(c, "".join(s)) for c, s in ref.items()])
for si, s in enumerate(SAMPLES):
    recs = [(f"{c}_h{h}", haplotype(c, si, h)) for c in CONTIGS for h in (0, 1)]
    write_fasta(f"{OUT}/{s}_hap.fa", recs)

with open(f"{OUT}/truth.tsv", "w") as fh:
    fh.write("chrom\tpos\tref\talt\t" + "\t".join(SAMPLES) + "\n")
    for c, pos, ra, alt, gts in sites:
        fh.write(f"{c}\t{pos + 1}\t{ra}\t{alt}\t" + "\t".join(map(str, gts)) + "\n")

with open(f"{OUT}/TruSeq3-PE.fa", "w") as fh:
    fh.write(">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n"
             ">PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n")

print(f"Simulated {len(sites)} variant sites across {len(SAMPLES)} samples -> {OUT}")
