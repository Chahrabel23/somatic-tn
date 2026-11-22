#!/usr/bin/env python3
"""Turn a VEP-annotated Mutect2 VCF into a flat table for review.

One line per variant. The transcript consequence kept is the one VEP flagged
with PICK=1 (falls back to the first CSQ entry). Tumour and normal sample names
are read from the ##tumor_sample / ##normal_sample header lines written by Mutect2.
"""
import argparse
import gzip


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_format(fmt, sample_field):
    return dict(zip(fmt.split(":"), sample_field.split(":")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--min-vaf", type=float, default=0.02)
    a = ap.parse_args()

    csq_fields, tumour, normal, samples = [], None, None, []
    rows = []
    with open_any(a.vcf) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = line.split("Format: ")[1].rstrip('">\n').split("|")
            elif line.startswith("##tumor_sample="):
                tumour = line.strip().split("=", 1)[1]
            elif line.startswith("##normal_sample="):
                normal = line.strip().split("=", 1)[1]
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                f = line.rstrip("\n").split("\t")
                info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
                csq = {}
                if "CSQ" in info and csq_fields:
                    entries = [dict(zip(csq_fields, e.split("|"))) for e in info["CSQ"].split(",")]
                    picked = [e for e in entries if e.get("PICK") == "1"]
                    csq = (picked or entries)[0]
                fmt = f[8]
                by_sample = {s: parse_format(fmt, f[9 + i]) for i, s in enumerate(samples)}
                t = by_sample.get(tumour, {})
                n = by_sample.get(normal, {})

                def ad_dp(d):
                    ad = [int(x) for x in d.get("AD", "0,0").split(",") if x not in (".", "")]
                    dp = sum(ad) if ad else 0
                    alt = ad[1] if len(ad) > 1 else 0
                    return alt, dp

                t_alt, t_dp = ad_dp(t)
                n_alt, n_dp = ad_dp(n)
                t_vaf = t_alt / t_dp if t_dp else 0.0
                n_vaf = n_alt / n_dp if n_dp else 0.0
                flag = "OK"
                if t_dp < a.min_depth:
                    flag = "LOW_DEPTH"
                elif t_vaf < a.min_vaf:
                    flag = "LOW_VAF"
                rows.append(
                    [
                        f[0], f[1], f[3], f[4],
                        csq.get("VARIANT_CLASS", ""),
                        csq.get("SYMBOL", ""),
                        csq.get("Feature", ""),
                        csq.get("Consequence", ""),
                        csq.get("IMPACT", ""),
                        csq.get("EXON", ""),
                        csq.get("HGVSc", ""),
                        csq.get("HGVSp", ""),
                        str(t_dp), str(t_alt), f"{t_vaf:.4f}",
                        str(n_dp), str(n_alt), f"{n_vaf:.4f}",
                        f[6], flag,
                    ]
                )

    header = [
        "chrom", "pos", "ref", "alt", "class", "gene", "transcript", "consequence",
        "impact", "exon", "hgvsc", "hgvsp", "tumour_depth", "tumour_alt", "tumour_vaf",
        "normal_depth", "normal_alt", "normal_vaf", "filter", "review_flag",
    ]
    with open(a.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"tumour={tumour} normal={normal} variants={len(rows)}")


if __name__ == "__main__":
    main()
