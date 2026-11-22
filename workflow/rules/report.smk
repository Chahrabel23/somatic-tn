# Per-pair HTML report and optional benchmark against a truth set


rule benchmark:
    input:
        calls="results/somatic/{pair}/{pair}.pass.vcf.gz",
        truth=lambda wc: TRUTH[wc.pair],
        fa="results/ref/genome.fa",
    output:
        tsv="results/benchmark/{pair}.benchmark.tsv",
        detail="results/benchmark/{pair}.benchmark_detail.tsv",
    log:
        "logs/benchmark/{pair}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools norm -m -both -f {input.fa} {input.truth} -Oz -o results/benchmark/{wildcards.pair}.truth.norm.vcf.gz 2> {log}
        python workflow/scripts/benchmark.py --calls {input.calls} \
          --truth results/benchmark/{wildcards.pair}.truth.norm.vcf.gz \
          --out {output.tsv} --detail {output.detail} >> {log} 2>&1
        """


rule cnv_benchmark:
    input:
        cns="results/cnv/{pair}/{pair}.call.cns",
        truth=lambda wc: TRUTH_CNV[wc.pair],
    output:
        "results/benchmark/{pair}.cnv_benchmark.tsv",
    log:
        "logs/benchmark/{pair}.cnv.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/cnv_benchmark.py --cns {input.cns} --truth {input.truth} \
          --out {output} > {log} 2>&1
        """


rule report:
    input:
        unpack(report_inputs),
    output:
        "results/report/{pair}.report.html",
    params:
        pair="{pair}",
        tumour=tumour_of,
        normal=normal_of,
        run_id=config.get("run_id", "RUN"),
        version=PIPELINE_VERSION,
        qc=config["qc"],
    log:
        "logs/report/{pair}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/report.py"
