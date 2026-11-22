# Each tool version is read inside the conda environment that runs the tool,
# so the report lists exactly what produced the results.


rule versions_align:
    output:
        "results/pipeline_info/versions/align.tsv",
    log:
        "logs/versions/align.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (printf "fastp\\t%s\\n" "$(fastp --version 2>&1 | cut -d' ' -f2)"
         printf "bwa\\t%s\\n" "$(bwa 2>&1 | grep Version | cut -d' ' -f2)"
         printf "samtools\\t%s\\n" "$(samtools --version | head -1 | cut -d' ' -f2)") > {output} 2> {log}
        """


rule versions_gatk:
    output:
        "results/pipeline_info/versions/gatk.tsv",
    log:
        "logs/versions/gatk.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        printf "gatk\\t%s\\n" "$(gatk --version 2>/dev/null | grep 'The Genome Analysis Toolkit' | rev | cut -d' ' -f1 | rev)" > {output} 2> {log}
        """


rule versions_bcftools:
    output:
        "results/pipeline_info/versions/bcftools.tsv",
    log:
        "logs/versions/bcftools.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        printf "bcftools\\t%s\\n" "$(bcftools --version | head -1 | cut -d' ' -f2)" > {output} 2> {log}
        """


rule versions_qc:
    output:
        "results/pipeline_info/versions/qc.tsv",
    log:
        "logs/versions/qc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        """
        (printf "mosdepth\\t%s\\n" "$(mosdepth --version | cut -d' ' -f2)"
         printf "multiqc\\t%s\\n" "$(multiqc --version | rev | cut -d' ' -f1 | rev)") > {output} 2> {log}
        """


rule versions_cnvkit:
    output:
        "results/pipeline_info/versions/cnvkit.tsv",
    log:
        "logs/versions/cnvkit.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        printf "cnvkit\\t%s\\n" "$(cnvkit.py version)" > {output} 2> {log}
        """


rule versions_vep:
    output:
        "results/pipeline_info/versions/vep.tsv",
    log:
        "logs/versions/vep.log",
    conda:
        "../envs/vep.yaml"
    shell:
        """
        printf "vep\\t%s\\n" "$(vep --help 2>&1 | grep 'ensembl-vep ' | tr -s ' ' | cut -d' ' -f4)" > {output} 2> {log}
        """
