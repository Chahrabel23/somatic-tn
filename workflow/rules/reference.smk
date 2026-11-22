# Reference preparation. The user's reference folder is never modified:
# everything is linked or written under results/ref/.


rule ref_fasta:
    input:
        REF["fasta"],
    output:
        fa="results/ref/genome.fa",
        fai="results/ref/genome.fa.fai",
        dict="results/ref/genome.dict",
    log:
        "logs/ref/fasta.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        ln -sf {input} {output.fa}
        samtools faidx {output.fa} 2> {log}
        gatk CreateSequenceDictionary -R {output.fa} -O {output.dict} >> {log} 2>&1
        """


rule bwa_index:
    input:
        "results/ref/genome.fa",
    output:
        multiext("results/ref/genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/ref/bwa_index.log",
    conda:
        "../envs/align.yaml"
    shell:
        "bwa index {input} > {log} 2>&1"


rule targets_bed:
    input:
        bed=REF["targets"],
        fai="results/ref/genome.fa.fai",
    output:
        "results/ref/targets.bed",
    log:
        "logs/ref/targets.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        # keep the first four columns, drop headers, sort in reference order
        """
        grep -vE '^(track|browser|#)' {input.bed} \
          | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,(NF>=4?$4:"-")}}' \
          | sort -k1,1 -k2,2n > {output} 2> {log}
        """


rule vep_gtf:
    input:
        REF["gtf"],
    output:
        gtf="results/ref/genes.gtf.gz",
        tbi="results/ref/genes.gtf.gz.tbi",
    log:
        "logs/ref/vep_gtf.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        (case {input} in *.gz) zcat {input};; *) cat {input};; esac) \
          | grep -v '^#' | sort -t$'\\t' -k1,1 -k4,4n -k5,5n \
          | bgzip -c > {output.gtf} 2> {log}
        tabix -p gff {output.gtf} >> {log} 2>&1
        """


rule cnvkit_access:
    input:
        "results/ref/genome.fa",
    output:
        "results/ref/access.bed",
    log:
        "logs/ref/cnvkit_access.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py access {input} -o {output} > {log} 2>&1"
