# Copy-number analysis with CNVkit, tumour against its matched normal.
# With a pool of normals, build reference.cnn once and pass it instead.


rule cnvkit_batch:
    input:
        t=lambda wc: recal_bam(tumour_of(wc)),
        n=lambda wc: recal_bam(normal_of(wc)),
        targets="results/ref/targets.bed",
        access="results/ref/access.bed",
        fa="results/ref/genome.fa",
    output:
        cnr="results/cnv/{pair}/{pair}.cnr",
        cns="results/cnv/{pair}/{pair}.cns",
        call="results/cnv/{pair}/{pair}.call.cns",
        ref="results/cnv/{pair}/{pair}.reference.cnn",
    params:
        outdir="results/cnv/{pair}/batch",
        method=PARAMS.get("cnvkit_method", "hybrid"),
        segment=PARAMS.get("cnvkit_segment", "cbs"),
        stem=tumour_of,  # CNVkit names files after the BAM stem (sample name)
    threads: 2
    log:
        "logs/cnvkit/{pair}.batch.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        cnvkit.py batch {input.t} --normal {input.n} --method {params.method} \
          --targets {input.targets} --access {input.access} --fasta {input.fa} \
          --segment-method {params.segment} --drop-low-coverage \
          --output-reference {output.ref} --output-dir {params.outdir} \
          --processes {threads} > {log} 2>&1
        cp {params.outdir}/{params.stem}.cnr {output.cnr}
        cp {params.outdir}/{params.stem}.cns {output.cns}
        cnvkit.py call {output.cns} -o {output.call} >> {log} 2>&1
        """


rule cnvkit_plots:
    input:
        cnr="results/cnv/{pair}/{pair}.cnr",
        cns="results/cnv/{pair}/{pair}.call.cns",
    output:
        scatter="results/cnv/{pair}/{pair}.scatter.png",
        genes="results/cnv/{pair}/{pair}.genemetrics.tsv",
    log:
        "logs/cnvkit/{pair}.plots.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        cnvkit.py scatter {input.cnr} -s {input.cns} -o {output.scatter} > {log} 2>&1
        cnvkit.py genemetrics {input.cnr} -s {input.cns} -t 0.2 -m 3 -o {output.genes} >> {log} 2>&1
        """


rule cnv_table:
    input:
        "results/cnv/{pair}/{pair}.call.cns",
    output:
        "results/cnv/{pair}/{pair}.cnv_segments.tsv",
    log:
        "logs/cnvkit/{pair}.table.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} NR==1{{for(i=1;i<=NF;i++) c[$i]=i;
             print "chromosome","start","end","genes","log2","copy_number","probes"; next}}
             {{print $c["chromosome"],$c["start"],$c["end"],$c["gene"],$c["log2"],$c["cn"],$c["probes"]}}' \
             {input} > {output} 2> {log}
        """
