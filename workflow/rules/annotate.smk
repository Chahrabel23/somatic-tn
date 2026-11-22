# Offline annotation with Ensembl VEP, using the GTF and FASTA directly
# (no VEP cache needed, so the run stays reproducible and offline).


rule vep:
    input:
        vcf="results/somatic/{pair}/{pair}.pass.vcf.gz",
        gtf="results/ref/genes.gtf.gz",
        tbi="results/ref/genes.gtf.gz.tbi",
        fa="results/ref/genome.fa",
        fai="results/ref/genome.fa.fai",
    output:
        vcf="results/annotation/{pair}/{pair}.vep.vcf.gz",
        summary="results/annotation/{pair}/{pair}.vep_summary.html",
    params:
        extra=PARAMS.get("vep", ""),
    threads: 2
    log:
        "logs/vep/{pair}.log",
    conda:
        "../envs/vep.yaml"
    shell:
        """
        vep --input_file {input.vcf} --output_file {output.vcf} --format vcf --vcf \
          --compress_output bgzip --gtf {input.gtf} --fasta {input.fa} \
          --stats_file {output.summary} --fork {threads} --force_overwrite \
          {params.extra} > {log} 2>&1
        """


rule variant_table:
    input:
        vcf="results/annotation/{pair}/{pair}.vep.vcf.gz",
    output:
        "results/somatic/{pair}/{pair}.variants.tsv",
    params:
        min_depth=config["qc"].get("min_tumour_depth", 100),
        min_vaf=config["qc"].get("min_tumour_vaf", 0.02),
    log:
        "logs/tables/{pair}.variants.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/variant_table.py --vcf {input.vcf} --out {output} \
          --min-depth {params.min_depth} --min-vaf {params.min_vaf} > {log} 2>&1
        """
