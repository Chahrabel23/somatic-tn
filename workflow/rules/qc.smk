# Coverage, alignment QC, MultiQC and traceability files


rule mosdepth:
    input:
        bam="results/bam/{sample}.recal.bam",
        bai="results/bam/{sample}.recal.bai",
        targets="results/ref/targets.bed",
    output:
        summary="results/qc/mosdepth/{sample}.mosdepth.summary.txt",
        dist="results/qc/mosdepth/{sample}.mosdepth.region.dist.txt",
        regions="results/qc/mosdepth/{sample}.regions.bed.gz",
        thresholds="results/qc/mosdepth/{sample}.thresholds.bed.gz",
    params:
        prefix="results/qc/mosdepth/{sample}",
        thresholds=PARAMS.get("mosdepth_thresholds", "100,250,500"),
    threads: 4
    log:
        "logs/mosdepth/{sample}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mosdepth --no-per-base --threads {threads} --by {input.targets} \
          --thresholds {params.thresholds} {params.prefix} {input.bam} > {log} 2>&1
        """


rule flagstat:
    input:
        "results/bam/{sample}.recal.bam",
    output:
        "results/qc/samtools/{sample}.flagstat.txt",
    log:
        "logs/samtools/{sample}.flagstat.log",
    conda:
        "../envs/align.yaml"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"


rule bcftools_stats:
    input:
        "results/somatic/{pair}/{pair}.pass.vcf.gz",
    output:
        "results/qc/bcftools/{pair}.bcftools_stats.txt",
    log:
        "logs/bcftools/{pair}.stats.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools stats {input} > {output} 2> {log}"


rule multiqc:
    input:
        expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        expand("results/qc/markdup/{sample}.metrics.txt", sample=SAMPLES),
        expand("results/qc/mosdepth/{sample}.mosdepth.summary.txt", sample=SAMPLES),
        expand("results/qc/samtools/{sample}.flagstat.txt", sample=SAMPLES),
        expand("results/qc/bcftools/{pair}.bcftools_stats.txt", pair=PAIRS),
        expand("results/annotation/{pair}/{pair}.vep_summary.html", pair=PAIRS),
    output:
        "results/multiqc/multiqc_report.html",
    params:
        title=config.get("run_id", "somatic-tn"),
    log:
        "logs/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc --force --title "{params.title}" -o results/multiqc \
          -n multiqc_report.html results/qc results/annotation > {log} 2>&1
        """


rule software_versions:
    input:
        expand("results/pipeline_info/versions/{env}.tsv", env=VERSION_ENVS),
    output:
        "results/pipeline_info/software_versions.tsv",
    params:
        version=PIPELINE_VERSION,
        smk=SNAKEMAKE_VERSION,
    log:
        "logs/versions/all.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (printf "tool\\tversion\\n"; printf "somatic-tn\\t{params.version}\\n";
         printf "snakemake\\t{params.smk}\\n"; cat {input}) > {output} 2> {log}
        """


rule input_checksums:
    input:
        fastq=[SAMPLES[s][k] for s in SAMPLES for k in ("fq1", "fq2")],
        ref=[REF["fasta"], REF["gtf"], REF["targets"], REF["germline_resource"], REF["common_snps"]]
        + list(REF["known_sites"]),
    output:
        md5="results/pipeline_info/input_checksums.md5",
        cfg="results/pipeline_info/config_used.yaml",
    log:
        "logs/versions/checksums.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/checksums.py"
