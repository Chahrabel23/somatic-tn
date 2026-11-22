# Read QC, mapping, duplicate marking and base quality recalibration


rule fastp:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["fq1"],
        r2=lambda wc: SAMPLES[wc.sample]["fq2"],
    output:
        r1=temp("results/fastp/{sample}.R1.fastq.gz"),
        r2=temp("results/fastp/{sample}.R2.fastq.gz"),
        json="results/qc/fastp/{sample}.fastp.json",
        html="results/qc/fastp/{sample}.fastp.html",
    params:
        extra=PARAMS.get("fastp", ""),
    threads: 4
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
          --json {output.json} --html {output.html} --report_title {wildcards.sample} \
          --thread {threads} {params.extra} > {log} 2>&1
        """


rule bwa_mem:
    input:
        r1="results/fastp/{sample}.R1.fastq.gz",
        r2="results/fastp/{sample}.R2.fastq.gz",
        idx=multiext("results/ref/genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        fa="results/ref/genome.fa",
    output:
        temp("results/bam/{sample}.sorted.bam"),
    params:
        rg=lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tLB:{wc.sample}\\tPL:ILLUMINA\\tPU:{wc.sample}",
        extra=PARAMS.get("bwa", ""),
    threads: 8
    log:
        "logs/bwa/{sample}.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (bwa mem -t {threads} {params.extra} -R '{params.rg}' {input.fa} {input.r1} {input.r2} \
          | samtools sort -@ 2 -m 1G -o {output} -) 2> {log}
        """


rule mark_duplicates:
    input:
        "results/bam/{sample}.sorted.bam",
    output:
        bam=temp("results/bam/{sample}.md.bam"),
        bai=temp("results/bam/{sample}.md.bai"),
        metrics="results/qc/markdup/{sample}.metrics.txt",
    params:
        java=JAVA,
    log:
        "logs/markdup/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" MarkDuplicates -I {input} -O {output.bam} \
          -M {output.metrics} --CREATE_INDEX true > {log} 2>&1
        """


rule base_recalibrator:
    input:
        bam="results/bam/{sample}.md.bam",
        bai="results/bam/{sample}.md.bai",
        fa="results/ref/genome.fa",
        dict="results/ref/genome.dict",
        known=REF["known_sites"],
        targets="results/ref/targets.bed",
    output:
        "results/bam/{sample}.recal.table",
    params:
        java=JAVA,
        known=known_sites_args(),
    log:
        "logs/bqsr/{sample}.recal.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" BaseRecalibrator -I {input.bam} -R {input.fa} \
          {params.known} -L {input.targets} --interval-padding 100 -O {output} > {log} 2>&1
        """


rule apply_bqsr:
    input:
        bam="results/bam/{sample}.md.bam",
        bai="results/bam/{sample}.md.bai",
        table="results/bam/{sample}.recal.table",
        fa="results/ref/genome.fa",
    output:
        bam="results/bam/{sample}.recal.bam",
        bai="results/bam/{sample}.recal.bai",
    params:
        java=JAVA,
    log:
        "logs/bqsr/{sample}.apply.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" ApplyBQSR -I {input.bam} -R {input.fa} \
          --bqsr-recal-file {input.table} -O {output.bam} --create-output-bam-index true > {log} 2>&1
        """
