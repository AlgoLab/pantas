rule get_transcripts:
    input:
        fa=FA,
        gtf=GTF,
    output:
        fa=pjoin(ODIR, "transcripts.fa"),
    conda:
        "../envs/gffread.yml"
    shell:
        """
        gffread -w {output.fa} -g {input.fa} {input.gtf}
        """


rule salmon_index:
    input:
        fa=rules.get_transcripts.output.fa,
    output:
        index=directory(pjoin(ODIR, "salmon", "salmon-index")),
    threads: workflow.cores
    conda:
        "../envs/salmon.yaml"
    log:
        pjoin(ODIR, "bench", "salmon", "salmon-index.time"),
    shell:
        """
        /usr/bin/time -vo {log} salmon index -t {input.fa} -i {output.index} -p {threads}
        """


rule salmon_quant:
    input:
        index=pjoin(ODIR, "salmon", "salmon-index"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
    output:
        outd=directory(pjoin(ODIR, "salmon", "{sample}/")),
        sf=pjoin(ODIR, "salmon", "{sample}/quant.sf"),
    threads: workflow.cores
    conda:
        "../envs/salmon.yaml"
    log:
        pjoin(ODIR, "bench", "salmon", "salmon-quant.{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log} salmon quant -i {input.index} -l IU -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output.outd} -p {threads}
        """


rule suppa2_generateevents:
    input:
        gtf=GTF,
    output:
        ioe=pjoin(ODIR, "suppa2", "annotated.events.ioe"),
    params:
        oprefix=pjoin(ODIR, "suppa2", "annotated"),
    threads: 1
    conda:
        "../envs/suppa2.yaml"
    log:
        pjoin(ODIR, "bench", "suppa2", "generateevents.time"),
    shell:
        """
        /usr/bin/time -vo {log} suppa.py generateEvents -i {input.gtf} -o {params.oprefix} -e SE SS RI -f ioe
        awk 'FNR==1 && NR!=1 {{ while (/^seqname/) getline; }} 1 {{print}}' {params.oprefix}*ioe > {output.ioe}
        """


rule suppa2:
    input:
        ioe=pjoin(ODIR, "suppa2", "annotated.events.ioe"),
        quants=expand(
            pjoin(ODIR, "salmon", "{sample}/quant.sf"),
            sample=FQs.keys(),
        ),
    output:
        outd=directory(pjoin(ODIR, "suppa2", "OUT")),
    params:
        quants=pjoin(ODIR, "salmon", "*/quant.sf"),
    threads: 1
    conda:
        "../envs/suppa2.yaml"
    log:
        pjoin(ODIR, "bench", "suppa2", "suppa2.time"),
    shell:
        """
        /usr/bin/time -vo {log} bash workflow/scripts/run_suppa2.sh {input.ioe} "{params.quants}" {output.outd}
        """
