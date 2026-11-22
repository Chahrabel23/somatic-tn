#!/usr/bin/env python3
"""Simulate a tumour-normal pair with known somatic variants (truth set).

Design
  * germline: 12 heterozygous SNVs present in BOTH normal and tumour (VAF ~50%).
    They must NOT be reported as somatic.
  * copy number: one gain (3 copies) and one loss (1 copy) in the tumour, written
    to truth_cnv.bed. VAFs of the somatic variants are kept the same inside them.
  * somatic: five groups of 13 variants (10 SNVs + 3 indels) at VAF 30%, 10%, 5%, 2% and 1%,
    present in the tumour only. Recall per group gives a first estimate of the limit of detection.

Reads are simulated with wgsim from haplotype FASTA files and mixed in the right
proportions. All variants fall inside the test target BED.

Usage: simulate_tumour_normal.py <ref.fa> <targets.bed> <outdir> [--seed 7]
"""
import argparse
import gzip
import random
import subprocess
from pathlib import Path

GENOME_LEN_FOR_DEPTH = None
READ_LEN, INSERT = 150, 350
VAFS = (0.30, 0.10, 0.05, 0.02, 0.01)
GAIN = (3500, 10500)    # inside the KCNMB3P1 target, tumour copy number 3
LOSS = (30200, 39600)   # inside the AP000547.3 target, tumour copy number 1


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq.append(line.strip().upper())
    return name, "".join(seq)


def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def apply(seq, variants):
    """variants: list of (pos1, ref, alt). Applied from the end so positions stay valid."""
    s = seq
    for pos, ref, alt in sorted(variants, key=lambda v: v[0], reverse=True):
        i = pos - 1
        assert s[i:i + len(ref)] == ref, (pos, ref, s[i:i + len(ref)])
        s = s[:i] + alt + s[i + len(ref):]
    return s


def wgsim(fa, n_pairs, prefix, seed):
    cmd = ["wgsim", "-N", str(n_pairs), "-1", str(READ_LEN), "-2", str(READ_LEN),
           "-d", str(INSERT), "-s", "30", "-e", "0.001", "-r", "0", "-R", "0", "-X", "0",
           "-S", str(seed), fa, f"{prefix}_1.fq", f"{prefix}_2.fq"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return f"{prefix}_1.fq", f"{prefix}_2.fq"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("targets")
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--tumour-depth", type=int, default=800)
    ap.add_argument("--normal-depth", type=int, default=400)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.outdir)
    (out / "tmp").mkdir(parents=True, exist_ok=True)

    chrom, seq = read_fasta(a.ref)
    targets = [(int(f[1]), int(f[2])) for f in (l.split("\t") for l in open(a.targets)) if f[0] == chrom]

    # candidate positions: inside targets, away from contig ends, no N, not in homopolymers
    def ok(p):
        i = p - 1
        if i < 400 or i > len(seq) - 400 or "N" in seq[i - 5:i + 10]:
            return False
        return not any(seq[j] == seq[j + 1] == seq[j + 2] for j in range(i - 3, i + 3))

    cand = [p for s, e in targets for p in range(s + 60, e - 60, 7) if ok(p)]
    rng.shuffle(cand)
    chosen = []
    for p in cand:
        if all(abs(p - q) > 200 for q in chosen):
            chosen.append(p)
        if len(chosen) == 12 + 13 * len(VAFS):
            break
    assert len(chosen) == 12 + 13 * len(VAFS), f"only {len(chosen)} positions found"

    def snv(p):
        ref = seq[p - 1]
        return (p, ref, rng.choice([b for b in "ACGT" if b != ref]))

    def indel(p, k):
        if k % 2 == 0:  # deletion of 1-6 bp
            n = rng.randint(1, 6)
            return (p, seq[p - 1:p + n], seq[p - 1])
        ins = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 5)))
        return (p, seq[p - 1], seq[p - 1] + ins)

    germ = [snv(p) for p in chosen[:12]]
    groups = {}
    rest = chosen[12:]
    for gi, vaf in enumerate(VAFS):
        ps = rest[gi * 13:(gi + 1) * 13]
        groups[vaf] = [snv(p) for p in ps[:10]] + [indel(p, k) for k, p in enumerate(ps[10:])]

    # haplotypes: H0 = reference, H1 = germline; each somatic group sits on H0 or H1 (half/half)
    variants_of = {"H0": [], "H1": list(germ)}
    for vaf, vs in groups.items():
        variants_of[f"H0_{vaf}"] = list(vs)
        variants_of[f"H1_{vaf}"] = germ + list(vs)

    def shift(pos, vs):
        return pos + sum(len(al) - len(r) for p, r, al in vs if p < pos)

    # tumour copy-number changes (normal stays diploid):
    #   GAIN: tandem duplication on the H0 haplotypes -> 3 copies
    #   LOSS: deletion on the H1 haplotypes          -> 1 copy (with LOH of germline SNVs)
    fa, length = {}, {}
    for k, vs in variants_of.items():
        s_hap = apply(seq, vs)
        for label, s_out in (("N", s_hap), ("T", s_hap)):
            if label == "T":
                if k.startswith("H0"):
                    g0, g1 = shift(GAIN[0], vs), shift(GAIN[1], vs)
                    s_out = s_hap[:g1] + s_hap[g0:g1] + s_hap[g1:]
                else:
                    l0, l1 = shift(LOSS[0], vs), shift(LOSS[1], vs)
                    s_out = s_hap[:l0] + s_hap[l1:]
            if label == "N" and k not in ("H0", "H1"):
                continue
            key = f"{label}_{k}"
            fa[key] = str(out / "tmp" / f"{key}.fa")
            length[key] = len(s_out)
            write_fasta(fa[key], chrom, s_out)

    pairs_per_x = len(seq) / (2 * READ_LEN)
    t_pairs = int(a.tumour_depth * pairs_per_x)
    n_pairs = int(a.normal_depth * pairs_per_x)

    mix_t = {}
    somatic_frac = sum(groups)  # 0.48
    for h in ("H0", "H1"):
        mix_t[h] = (1 - somatic_frac) / 2
        for vaf in groups:
            mix_t[f"{h}_{vaf}"] = vaf / 2
    mix_n = {"H0": 0.5, "H1": 0.5}
    mix_t = {f"T_{k}": v for k, v in mix_t.items()}
    mix_n = {f"N_{k}": v for k, v in mix_n.items()}

    def simulate(mix, total, label, seed0):
        r1s, r2s = [], []
        for i, (h, frac) in enumerate(sorted(mix.items())):
            n = max(1, int(total * frac * length[h] / len(seq)))  # keep depth per copy constant
            r1, r2 = wgsim(fa[h], n, str(out / "tmp" / f"{label}_{h}"), seed0 + i)
            r1s.append(r1)
            r2s.append(r2)
        # wgsim writes a flat Q17 for every base; the simulated error rate is 0.1%,
        # so the qualities are set to Q30 to match it.
        for reads, mate in ((r1s, 1), (r2s, 2)):
            target = out / f"{label}_R{mate}.fastq.gz"
            with gzip.open(target, "wt", compresslevel=4) as fh:
                for path in reads:
                    with open(path) as src:
                        for i, line in enumerate(src):
                            fh.write("?" * (len(line) - 1) + "\n" if i % 4 == 3 else line)

    simulate(mix_t, t_pairs, "SIM_T", a.seed * 100)
    simulate(mix_n, n_pairs, "SIM_N", a.seed * 100 + 50)

    # truth VCF (somatic only) and germline list
    def write_vcf(path, recs):
        with open(path, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write(f"##contig=<ID={chrom},length={len(seq)}>\n")
            fh.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Simulated variant allele fraction">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for (p, r, al), vaf in sorted(recs, key=lambda x: x[0][0]):
                fh.write(f"{chrom}\t{p}\t.\t{r}\t{al}\t.\tPASS\tVAF={vaf}\n")

    write_vcf(out / "truth_somatic.vcf", [(v, vaf) for vaf, vs in groups.items() for v in vs])
    write_vcf(out / "germline_het.vcf", [(v, 0.5) for v in germ])
    with open(out / "truth_cnv.bed", "w") as fh:
        fh.write(f"{chrom}\t{GAIN[0]}\t{GAIN[1]}\tGAIN\t3\n")
        fh.write(f"{chrom}\t{LOSS[0]}\t{LOSS[1]}\tLOSS\t1\n")
    for f in (out / "tmp").iterdir():
        f.unlink()
    (out / "tmp").rmdir()
    print(f"tumour pairs={t_pairs} normal pairs={n_pairs} somatic={sum(len(v) for v in groups.values())} germline={len(germ)}")


if __name__ == "__main__":
    main()
