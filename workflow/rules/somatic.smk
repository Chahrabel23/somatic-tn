# Somatic SNV / indel calling with GATK Mutect2 (tumour-normal mode)
# Mutect2 -> LearnReadOrientationModel -> GetPileupSummaries ->
# CalculateContamination -> FilterMutectCalls -> normalise + PASS


rule mutect2:
    input:
        t=lambda wc: recal_bam(tumour_of(wc)),
        n=lambda wc: recal_bam(normal_of(wc)),
        fa="results/ref/genome.fa",
        dict="results/ref/genome.dict",
        targets="results/ref/targets.bed",
        germline=REF["germline_resource"],
    output:
        vcf="results/somatic/{pair}/{pair}.mutect2.unfiltered.vcf.gz",
        stats="results/somatic/{pair}/{pair}.mutect2.unfiltered.vcf.gz.stats",
        f1r2="results/somatic/{pair}/{pair}.f1r2.tar.gz",
        bam="results/somatic/{pair}/{pair}.mutect2.bamout.bam",
    params:
        java=JAVA,
        normal=normal_of,
        pon=pon_arg(),
        extra=PARAMS.get("mutect2", ""),
    threads: 4
    log:
        "logs/mutect2/{pair}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" Mutect2 -R {input.fa} \
          -I {input.t} -I {input.n} -normal {params.normal} \
          --germline-resource {input.germline} {params.pon} \
          -L {input.targets} --interval-padding 50 \
          --native-pair-hmm-threads {threads} \
          --f1r2-tar-gz {output.f1r2} --bam-output {output.bam} \
          {params.extra} -O {output.vcf} > {log} 2>&1
        """


rule learn_read_orientation:
    input:
        "results/somatic/{pair}/{pair}.f1r2.tar.gz",
    output:
        "results/somatic/{pair}/{pair}.read-orientation-model.tar.gz",
    params:
        java=JAVA,
    log:
        "logs/mutect2/{pair}.orientation.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" LearnReadOrientationModel -I {input} -O {output} > {log} 2>&1
        """


rule pileup_summaries:
    input:
        bam="results/bam/{sample}.recal.bam",
        common=REF["common_snps"],
        targets="results/ref/targets.bed",
    output:
        "results/somatic/pileups/{sample}.pileups.table",
    params:
        java=JAVA,
    log:
        "logs/pileups/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" GetPileupSummaries -I {input.bam} \
          -V {input.common} -L {input.targets} -O {output} > {log} 2>&1
        """


rule calculate_contamination:
    input:
        t=lambda wc: f"results/somatic/pileups/{tumour_of(wc)}.pileups.table",
        n=lambda wc: f"results/somatic/pileups/{normal_of(wc)}.pileups.table",
    output:
        cont="results/somatic/{pair}/{pair}.contamination.table",
        seg="results/somatic/{pair}/{pair}.segments.table",
    params:
        java=JAVA,
    log:
        "logs/mutect2/{pair}.contamination.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" CalculateContamination -I {input.t} -matched {input.n} \
          -O {output.cont} --tumor-segmentation {output.seg} > {log} 2>&1
        """


rule filter_mutect:
    input:
        vcf="results/somatic/{pair}/{pair}.mutect2.unfiltered.vcf.gz",
        stats="results/somatic/{pair}/{pair}.mutect2.unfiltered.vcf.gz.stats",
        ob="results/somatic/{pair}/{pair}.read-orientation-model.tar.gz",
        cont="results/somatic/{pair}/{pair}.contamination.table",
        seg="results/somatic/{pair}/{pair}.segments.table",
        fa="results/ref/genome.fa",
    output:
        "results/somatic/{pair}/{pair}.mutect2.filtered.vcf.gz",
    params:
        java=JAVA,
        extra=PARAMS.get("filter_mutect", ""),
    log:
        "logs/mutect2/{pair}.filter.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk --java-options "{params.java}" FilterMutectCalls -R {input.fa} -V {input.vcf} \
          --stats {input.stats} --ob-priors {input.ob} \
          --contamination-table {input.cont} --tumor-segmentation {input.seg} \
          {params.extra} -O {output} > {log} 2>&1
        """


rule normalise_pass:
    input:
        vcf="results/somatic/{pair}/{pair}.mutect2.filtered.vcf.gz",
        fa="results/ref/genome.fa",
    output:
        vcf="results/somatic/{pair}/{pair}.pass.vcf.gz",
        tbi="results/somatic/{pair}/{pair}.pass.vcf.gz.tbi",
    log:
        "logs/bcftools/{pair}.norm.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        (bcftools norm -m -both -f {input.fa} {input.vcf} \
          | bcftools view -f PASS -Oz -o {output.vcf}) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """
