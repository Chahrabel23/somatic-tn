#!/usr/bin/env python3
"""Compare PASS somatic calls with a truth VCF.

Both files must be normalised (bcftools norm -m -both -f ref) beforehand.
A call is a true positive when CHROM, POS, REF and ALT match exactly.
If the truth VCF has an INFO/VAF tag, recall is also reported per VAF level,
which gives a first estimate of the limit of detection.
"""
import argparse
import gzip
from collections import defaultdict


def records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
            kind = "SNV" if len(f[3]) == 1 and len(f[4]) == 1 else "INDEL"
            yield (f[0], int(f[1]), f[3], f[4]), kind, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--detail", required=True)
    a = ap.parse_args()

    truth = {k: (kind, info.get("VAF", "NA")) for k, kind, info in records(a.truth)}
    calls = {k: kind for k, kind, _ in records(a.calls)}

    tp = {k for k in calls if k in truth}
    fp = {k for k in calls if k not in truth}
    fn = {k for k in truth if k not in calls}

    def ratio(x, y):
        return f"{x / y:.3f}" if y else "NA"

    lines = [["stratum", "truth", "called_tp", "fn", "fp", "recall", "precision"]]
    for kind in ("ALL", "SNV", "INDEL"):
        sel = (lambda k, kd: True) if kind == "ALL" else (lambda k, kd, kind=kind: kd == kind)
        t = [k for k, (kd, _) in truth.items() if sel(k, kd)]
        c_tp = [k for k in tp if sel(k, truth[k][0])]
        c_fp = [k for k in fp if sel(k, calls[k])]
        lines.append([kind, len(t), len(c_tp), len(t) - len(c_tp), len(c_fp),
                      ratio(len(c_tp), len(t)), ratio(len(c_tp), len(c_tp) + len(c_fp))])

    by_vaf = defaultdict(lambda: [0, 0])
    for k, (kd, vaf) in truth.items():
        by_vaf[vaf][0] += 1
        by_vaf[vaf][1] += k in tp
    for vaf in sorted(by_vaf, key=lambda v: float(v) if v != "NA" else -1, reverse=True):
        n, hit = by_vaf[vaf]
        lines.append([f"VAF={vaf}", n, hit, n - hit, "", ratio(hit, n), ""])

    with open(a.out, "w") as out:
        for row in lines:
            out.write("\t".join(map(str, row)) + "\n")
    with open(a.detail, "w") as out:
        out.write("chrom\tpos\tref\talt\ttype\ttruth_vaf\tstatus\n")
        for k in sorted(set(truth) | set(calls), key=lambda x: (x[0], x[1])):
            status = "TP" if k in tp else ("FP" if k in fp else "FN")
            kind = truth[k][0] if k in truth else calls[k]
            vaf = truth[k][1] if k in truth else ""
            out.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{kind}\t{vaf}\t{status}\n")
    print("\n".join("\t".join(map(str, r)) for r in lines))


if __name__ == "__main__":
    main()
