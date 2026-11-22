import csv
import os
from pathlib import Path


def _abs(path):
    return str(Path(path).expanduser().resolve()) if path else ""


def read_samples(path):
    with open(path) as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r.get("sample")]
    needed = {"patient", "sample", "status", "fq1", "fq2"}
    missing = needed - set(rows[0].keys()) if rows else needed
    if missing:
        raise ValueError(f"Sample sheet {path} is missing columns: {sorted(missing)}")
    table = {}
    for r in rows:
        status = r["status"].strip().lower().replace("tumor", "tumour")
        if status not in {"tumour", "normal"}:
            raise ValueError(f"Sample {r['sample']}: status must be tumour or normal, got {r['status']}")
        if r["sample"] in table:
            raise ValueError(f"Duplicated sample name: {r['sample']}")
        table[r["sample"]] = {
            "patient": r["patient"].strip(),
            "status": status,
            "fq1": _abs(r["fq1"].strip()),
            "fq2": _abs(r["fq2"].strip()),
        }
    return table


def make_pairs(samples):
    pairs = {}
    patients = sorted({s["patient"] for s in samples.values()})
    for p in patients:
        normals = [n for n, s in samples.items() if s["patient"] == p and s["status"] == "normal"]
        tumours = [t for t, s in samples.items() if s["patient"] == p and s["status"] == "tumour"]
        if len(normals) != 1:
            raise ValueError(f"Patient {p}: expected exactly one normal sample, found {len(normals)}")
        for t in tumours:
            pairs[f"{t}_vs_{normals[0]}"] = {"tumour": t, "normal": normals[0], "patient": p}
    if not pairs:
        raise ValueError("No tumour-normal pair could be built from the sample sheet")
    return pairs


SAMPLES = read_samples(config["samples"])
PAIRS = make_pairs(SAMPLES)
TRUTH = {k: _abs(v) for k, v in (config.get("truth") or {}).items() if k in PAIRS}
TRUTH_CNV = {k: _abs(v) for k, v in (config.get("truth_cnv") or {}).items() if k in PAIRS}

REF = {k: _abs(v) if isinstance(v, str) else [_abs(x) for x in v] for k, v in config["reference"].items()}
PARAMS = config["params"]
JAVA = PARAMS.get("java_opts", "-Xmx4g")

wildcard_constraints:
    sample="|".join(sorted(SAMPLES)),
    pair="|".join(sorted(PAIRS)),


def tumour_of(wc):
    return PAIRS[wc.pair]["tumour"]


def normal_of(wc):
    return PAIRS[wc.pair]["normal"]


def recal_bam(sample):
    return f"results/bam/{sample}.recal.bam"


def known_sites_args():
    return " ".join(f"--known-sites {k}" for k in REF["known_sites"])


def pon_arg():
    pon = REF.get("panel_of_normals") or ""
    return f"--panel-of-normals {pon}" if pon else ""


VERSION_ENVS = ["align", "gatk", "bcftools", "qc", "cnvkit", "vep"]


def report_inputs(wc):
    t, n = tumour_of(wc), normal_of(wc)
    files = {
        "variants": f"results/somatic/{wc.pair}/{wc.pair}.variants.tsv",
        "contamination": f"results/somatic/{wc.pair}/{wc.pair}.contamination.table",
        "filtered": f"results/somatic/{wc.pair}/{wc.pair}.mutect2.filtered.vcf.gz",
        "cnv": f"results/cnv/{wc.pair}/{wc.pair}.cnv_segments.tsv",
        "scatter": f"results/cnv/{wc.pair}/{wc.pair}.scatter.png",
        "versions": "results/pipeline_info/software_versions.tsv",
        "checksums": "results/pipeline_info/input_checksums.md5",
    }
    for role, s in (("tumour", t), ("normal", n)):
        files[f"{role}_fastp"] = f"results/qc/fastp/{s}.fastp.json"
        files[f"{role}_dup"] = f"results/qc/markdup/{s}.metrics.txt"
        files[f"{role}_mosdepth"] = f"results/qc/mosdepth/{s}.mosdepth.summary.txt"
        files[f"{role}_thresholds"] = f"results/qc/mosdepth/{s}.thresholds.bed.gz"
    if wc.pair in TRUTH:
        files["benchmark"] = f"results/benchmark/{wc.pair}.benchmark.tsv"
    if wc.pair in TRUTH_CNV:
        files["cnv_benchmark"] = f"results/benchmark/{wc.pair}.cnv_benchmark.tsv"
    return files
