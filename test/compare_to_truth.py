#!/usr/bin/env python3
"""Compare the pipeline's PASS variants with the simulated truth (site level + genotypes)."""
import gzip
import os

HERE = os.path.dirname(os.path.abspath(__file__))
truth_path = os.path.join(HERE, "data", "truth.tsv")
vcf_path = os.path.join(HERE, "results", "variants", "test_cohort.PASS.vcf.gz")

# Reference sequence, needed to left-align indels the same way GATK does
REF, name = {}, None
with open(os.path.join(HERE, "data", "ref.fa")) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            REF[name] = []
        else:
            REF[name].append(line)
REF = {k: "".join(v) for k, v in REF.items()}


def normalize(c, pos, ref, alt):
    """Left-align and trim a variant (Tan et al. 2015), so equivalent indels compare equal."""
    seq = REF[c]
    changed = True
    while changed:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt, changed = ref[:-1], alt[:-1], True
        if not ref or not alt:
            pos -= 1
            base = seq[pos - 1]
            ref, alt, changed = base + ref, base + alt, True
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return (c, pos, ref, alt)


truth = {}
with open(truth_path) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    samples = header[4:]
    for line in fh:
        c, pos, ref, alt, *gts = line.rstrip("\n").split("\t")
        truth[normalize(c, int(pos), ref, alt)] = [int(g) for g in gts]

called, vcf_samples = {}, []
with gzip.open(vcf_path, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            vcf_samples = f[9:]
            continue
        gt_idx = f[8].split(":").index("GT")
        gts = []
        for s in f[9:]:
            gt = s.split(":")[gt_idx].replace("|", "/")
            gts.append(None if "." in gt else sum(a != "0" for a in gt.split("/")))
        for alt in f[4].split(","):
            called[normalize(f[0], int(f[1]), f[3], alt)] = gts

kind = lambda k: "SNP" if len(k[2]) == 1 and len(k[3]) == 1 else "INDEL"
print(f"{'type':6} {'truth':>6} {'called':>7} {'found':>6} {'recall':>7} {'precision':>9}")
for t in ("SNP", "INDEL"):
    T = {k for k in truth if kind(k) == t}
    C = {k for k in called if kind(k) == t}
    tp = len(T & C)
    print(f"{t:6} {len(T):6d} {len(C):7d} {tp:6d} {tp / max(len(T), 1):7.1%} {tp / max(len(C), 1):9.1%}")

order = [vcf_samples.index(s) for s in samples]
agree = total = 0
for k, tg in truth.items():
    if k in called:
        cg = [called[k][i] for i in order]
        for a, b in zip(tg, cg):
            if b is not None:
                total += 1
                agree += (a == b)
print(f"Genotype concordance at shared sites: {agree}/{total} = {agree / max(total, 1):.1%}")
