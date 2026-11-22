"""Build a self-contained HTML report for one tumour-normal pair.

Called by Snakemake (rule report). Uses only the Python standard library.
"""
import base64
import csv
import datetime
import gzip
import html
import json
from collections import Counter

sm = snakemake  # noqa: F821  (injected by Snakemake)
inp, par = sm.input, sm.params
QC = par.qc


def esc(x):
    return html.escape(str(x))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fastp_stats(path):
    d = json.load(open(path))["summary"]
    return {
        "reads_raw": d["before_filtering"]["total_reads"],
        "reads_kept": d["after_filtering"]["total_reads"],
        "q30": d["after_filtering"]["q30_rate"],
        "read_len": d["after_filtering"]["read1_mean_length"],
    }


def dup_rate(path):
    lines = [l for l in open(path) if l.strip() and not l.startswith("#")]
    head = lines[0].rstrip("\n").split("\t")
    vals = lines[1].rstrip("\n").split("\t")
    return float(dict(zip(head, vals))["PERCENT_DUPLICATION"] or 0)


def mean_cov(path):
    for row in read_tsv(path):
        if row["chrom"] == "total_region":
            return float(row["mean"])
    rows = read_tsv(path)
    return float(rows[-1]["mean"]) if rows else 0.0


def pct_above(path):
    total, sums, labels = 0, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                labels = f[4:]
                sums = [0] * len(labels)
                continue
            total += int(f[2]) - int(f[1])
            for i, v in enumerate(f[4:]):
                sums[i] += int(v)
    return {lab.lower(): (100 * s / total if total else 0) for lab, s in zip(labels, sums)}


def change(ref, alt):
    if len(ref) <= 10 and len(alt) <= 10:
        return f"{ref}&gt;{alt}"
    d = len(alt) - len(ref)
    return f"{'ins' if d > 0 else 'del'} {abs(d)} bp"


def badge(ok):
    return '<span class="ok">PASS</span>' if ok else '<span class="warn">CHECK</span>'


# ---------- sample QC ----------
qc_rows = []
for role, name in (("tumour", par.tumour), ("normal", par.normal)):
    fp = fastp_stats(inp[f"{role}_fastp"])
    dup = dup_rate(inp[f"{role}_dup"])
    cov = mean_cov(inp[f"{role}_mosdepth"])
    above = pct_above(inp[f"{role}_thresholds"])
    pct250 = above.get("250x", 0.0)
    checks = [
        cov >= QC["min_mean_coverage"],
        pct250 >= QC["min_pct_target_250x"],
        dup <= QC["max_duplication"],
        fp["q30"] >= QC["min_q30_rate"],
    ]
    qc_rows.append((role, name, fp, dup, cov, above, all(checks)))

# ---------- contamination ----------
cont = read_tsv(inp.contamination)[0]
cont_val = float(cont["contamination"])
cont_err = float(cont["error"])

# ---------- Mutect2 filter summary ----------
filters = Counter()
with gzip.open(inp.filtered, "rt") as fh:
    for line in fh:
        if not line.startswith("#"):
            for flt in line.split("\t")[6].split(";"):
                filters[flt] += 1

variants = read_tsv(inp.variants)
cnv = read_tsv(inp.cnv)
versions = read_tsv(inp.versions)
checksums = [l.split(None, 1) for l in open(inp.checksums) if l.strip()]
bench = read_tsv(inp.benchmark) if "benchmark" in inp.keys() else None
cnv_bench = read_tsv(inp.cnv_benchmark) if "cnv_benchmark" in inp.keys() else None
png = base64.b64encode(open(inp.scatter, "rb").read()).decode()

# ---------- HTML ----------
css = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;color:#1d2433;background:#f6f7fb}
header{background:#3b2a63;color:#fff;padding:22px 36px}
header h1{margin:0;font-size:22px} header p{margin:4px 0 0;opacity:.85;font-size:13px}
main{padding:24px 36px;max-width:1200px}
section{background:#fff;border:1px solid #e2e5ee;border-radius:8px;padding:16px 20px;margin-bottom:18px}
h2{font-size:16px;margin:0 0 10px;color:#3b2a63}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{border-bottom:1px solid #eceef4;padding:6px 8px;text-align:left;vertical-align:top}
th{background:#f1f0f7;font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.ok{background:#e3f4e8;color:#166534;padding:2px 7px;border-radius:10px;font-weight:600;font-size:11px}
.warn{background:#fdf1dc;color:#92400e;padding:2px 7px;border-radius:10px;font-weight:600;font-size:11px}
.muted{color:#6b7280;font-size:12px}
img{max-width:100%;border:1px solid #eceef4;border-radius:6px}
code{background:#f1f0f7;padding:1px 4px;border-radius:3px;font-size:11.5px}
"""

out = []
w = out.append
w(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{esc(par.pair)} somatic report</title>")
w(f"<style>{css}</style></head><body>")
w(f"<header><h1>Somatic tumour-normal report: {esc(par.pair)}</h1>")
w(f"<p>Run {esc(par.run_id)} · tumour {esc(par.tumour)} · normal {esc(par.normal)} · "
  f"somatic-tn v{esc(par.version)} · generated {datetime.datetime.now():%Y-%m-%d %H:%M}</p></header><main>")

# QC
w("<section><h2>1. Sample QC</h2><table><tr><th>Role</th><th>Sample</th><th>Reads (raw)</th>"
  "<th>Reads kept</th><th>Q30</th><th>Duplication</th><th>Mean target depth</th>")
labels = list(qc_rows[0][5].keys())
for lab in labels:
    w(f"<th>Target ≥{esc(lab)}</th>")
w("<th>Status</th></tr>")
for role, name, fp, dup, cov, above, ok in qc_rows:
    w(f"<tr><td>{role}</td><td>{esc(name)}</td><td class='num'>{fp['reads_raw']:,}</td>"
      f"<td class='num'>{fp['reads_kept']:,}</td><td class='num'>{fp['q30']*100:.1f}%</td>"
      f"<td class='num'>{dup*100:.1f}%</td><td class='num'>{cov:,.0f}x</td>")
    for lab in labels:
        w(f"<td class='num'>{above[lab]:.1f}%</td>")
    w(f"<td>{badge(ok)}</td></tr>")
w("</table>")
w(f"<p class='muted'>Criteria: mean depth ≥ {QC['min_mean_coverage']}x, ≥ {QC['min_pct_target_250x']}% of target bases "
  f"at 250x, duplication ≤ {QC['max_duplication']*100:.0f}%, Q30 ≥ {QC['min_q30_rate']*100:.0f}%. "
  "Thresholds are set in the config file and must come from the lab's own validation.</p></section>")

# contamination + filters
w("<section><h2>2. Contamination and Mutect2 filtering</h2>")
w(f"<p>Estimated cross-sample contamination in the tumour: <b>{cont_val*100:.2f}%</b> "
  f"(± {cont_err*100:.2f}%) {badge(cont_val <= QC['max_contamination'])}</p>")
w("<table><tr><th>FILTER</th><th>Variants</th></tr>")
for flt, n in filters.most_common():
    w(f"<tr><td><code>{esc(flt)}</code></td><td class='num'>{n}</td></tr>")
w("</table></section>")

# variants
w(f"<section><h2>3. Somatic SNVs and indels (PASS, n = {len(variants)})</h2>")
if variants:
    w("<table><tr><th>Position</th><th>Change</th><th>Gene</th><th>Consequence</th><th>HGVSc</th>"
      "<th>HGVSp</th><th>Tumour VAF</th><th>Tumour depth</th><th>Normal VAF</th><th>Review</th></tr>")
    for v in variants:
        flag = '<span class="ok">OK</span>' if v["review_flag"] == "OK" else f'<span class="warn">{esc(v["review_flag"])}</span>'
        w(f"<tr><td>{esc(v['chrom'])}:{esc(v['pos'])}</td><td>{change(esc(v['ref']), esc(v['alt']))}</td>"
          f"<td>{esc(v['gene'])}</td><td>{esc(v['consequence'])}</td><td>{esc(v['hgvsc'].split(':')[-1])}</td>"
          f"<td>{esc(v['hgvsp'].split(':')[-1])}</td><td class='num'>{float(v['tumour_vaf'])*100:.1f}%</td>"
          f"<td class='num'>{esc(v['tumour_depth'])}</td><td class='num'>{float(v['normal_vaf'])*100:.1f}%</td>"
          f"<td>{flag}</td></tr>")
    w("</table>")
else:
    w("<p>No variant passed the filters.</p>")
w("<p class='muted'>Annotation: Ensembl VEP offline (GTF + FASTA). Clinical classification "
  "(e.g. AMP/ASCO/CAP tiers) is not done by the pipeline and stays with the reviewing scientist.</p></section>")

# CNV
w("<section><h2>4. Copy number (CNVkit)</h2>")
w(f"<img src='data:image/png;base64,{png}' alt='CNVkit scatter plot'>")
w("<table><tr><th>Segment</th><th>Genes</th><th>log2</th><th>Copy number</th><th>Bins</th></tr>")
for s in cnv:
    genes = s["genes"] if len(s["genes"]) < 80 else s["genes"][:77] + "..."
    w(f"<tr><td>{esc(s['chromosome'])}:{esc(s['start'])}-{esc(s['end'])}</td><td>{esc(genes)}</td>"
      f"<td class='num'>{float(s['log2']):.2f}</td><td class='num'>{esc(s['copy_number'])}</td>"
      f"<td class='num'>{esc(s['probes'])}</td></tr>")
w("</table><p class='muted'>Copy number assumes a diploid genome and 100% purity unless purity "
  "is given to <code>cnvkit.py call</code>.</p></section>")

# benchmark
def table(rows):
    w("<table><tr>" + "".join(f"<th>{esc(k)}</th>" for k in rows[0].keys()) + "</tr>")
    for r in rows:
        w("<tr>" + "".join(f"<td>{esc(v)}</td>" for v in r.values()) + "</tr>")
    w("</table>")


if bench or cnv_bench:
    w("<section><h2>5. Benchmark against the truth set</h2>")
    if bench:
        w("<p>SNVs and indels (exact match on position, REF and ALT after normalisation)</p>")
        table(bench)
    if cnv_bench:
        w("<p>Copy-number regions</p>")
        table(cnv_bench)
    w("</section>")

# traceability
w("<section><h2>Traceability</h2><table><tr><th>Tool</th><th>Version</th></tr>")
for v in versions:
    w(f"<tr><td>{esc(v['tool'])}</td><td>{esc(v['version'])}</td></tr>")
w("</table><p class='muted'>Input MD5 checksums</p><table><tr><th>MD5</th><th>File</th></tr>")
for md5, path in checksums:
    w(f"<tr><td><code>{esc(md5)}</code></td><td>{esc(path.strip().split('/')[-1])}</td></tr>")
w("</table></section>")

w("<section><h2>Scope</h2><p class='muted'>Research use only. This workflow has not been validated "
  "for diagnostic use. Gene fusions are not called from DNA here; an RNA module (Arriba / STAR-Fusion) "
  "is planned.</p></section>")
w("</main></body></html>")

with open(sm.output[0], "w") as fh:
    fh.write("\n".join(out))
