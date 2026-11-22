#!/usr/bin/env python3
"""Compare CNVkit called segments (.call.cns) with a truth BED.

Truth BED columns: chrom, start, end, label, expected copy number.
For each truth region the segment covering most of it is taken. Called segments
with copy number != 2 that do not overlap any truth region are listed as
unexpected calls.
"""
import argparse
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cns", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-fraction", type=float, default=0.5,
                    help="minimum fraction of the truth region covered by the matching segment")
    a = ap.parse_args()

    segs = list(csv.DictReader(open(a.cns), delimiter="\t"))
    truth = [l.rstrip("\n").split("\t") for l in open(a.truth) if l.strip()]

    rows = [["region", "label", "expected_cn", "called_cn", "log2", "overlap_fraction", "status"]]
    used = set()
    for chrom, start, end, label, cn in truth:
        start, end = int(start), int(end)
        best, best_ov = None, 0
        for i, s in enumerate(segs):
            if s["chromosome"] != chrom:
                continue
            ov = min(end, int(s["end"])) - max(start, int(s["start"]))
            if ov > best_ov:
                best, best_ov = i, ov
        if best is None:
            rows.append([f"{chrom}:{start}-{end}", label, cn, "NA", "NA", "0", "MISSED"])
            continue
        used.add(best)
        s = segs[best]
        frac = best_ov / (end - start)
        status = "MATCH" if s["cn"] == cn and frac >= a.min_fraction else "MISMATCH"
        rows.append([f"{chrom}:{start}-{end}", label, cn, s["cn"], f"{float(s['log2']):.2f}",
                     f"{frac:.2f}", status])
    for i, s in enumerate(segs):
        if i not in used and s["cn"] != "2":
            rows.append([f"{s['chromosome']}:{s['start']}-{s['end']}", "unexpected", "2", s["cn"],
                         f"{float(s['log2']):.2f}", "", "FALSE_CALL"])
    with open(a.out, "w") as out:
        for r in rows:
            out.write("\t".join(r) + "\n")
    print("\n".join("\t".join(r) for r in rows))


if __name__ == "__main__":
    main()
