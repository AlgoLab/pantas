rule get_transcripts:
    input:
        fa=FA,
        gtf=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf"),
    output:
        fa=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.cdna.fa"),
    conda:
        pjoin(ENVS, "gffread.yaml")
    shell:
        """
        gffread -w {output.fa} -g {input.fa} {input.gtf}
        """


rule salmon_index:
    input:
        fa=rules.get_transcripts.output.fa,
    output:
        index=directory(pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-index")),
    threads: workflow.cores
    conda:
        pjoin(ENVS, "salmon.yaml")
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "salmon", "index.time"),
    shell:
        """
        /usr/bin/time -vo {log} salmon index -t {input.fa} -i {output.index} -p {threads}
        """


rule salmon_quant:
    input:
        index=pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-index"),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    output:
        outd=directory(pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-quant-{x}/")),
        sf=pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-quant-{x}", "quant.sf"),
    threads: workflow.cores
    conda:
        pjoin(ENVS, "salmon.yaml")
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "salmon", "quant-{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} salmon quant -i {input.index} -l IU -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output.outd} -p {threads}
        """
