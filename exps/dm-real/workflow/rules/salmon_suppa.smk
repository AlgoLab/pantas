rule salmon_index:
    input:
        fa=FA
    output:
        index=directory(pjoin(ODIR, "salmon_suppa", "salmon-index.k{k}")),
    threads: workflow.cores
    conda: "../envs/salmon.yaml"
    log: pjoin(ODIR , "bench", "salmon_suppa", "salmon-index.k{k}.time")
    shell:
        """
        /usr/bin/time -vo {log} salmon index -t {input.fa} -i {output.index} -p {threads} -k{wildcards.k}
        """

rule salmon_quant:
    input:
        index=pjoin(ODIR, "salmon_suppa", "salmon-index.k{k}"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
    output:
        outd=directory(pjoin(ODIR, "salmon_suppa", "{sample}.salmon.k{k}/")),
        sf=pjoin(ODIR, "salmon_suppa", "{sample}.salmon.k{k}/quant.sf"),
    threads: workflow.cores
    conda: "../envs/salmon.yaml"
    log: pjoin(ODIR, "bench", "salmon_suppa", "salmon-quant.k{k}.{sample}.time")
    shell:
        """
        /usr/bin/time -vo {log} salmon quant -i {input.index} -l IU -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output.outd} -p {threads}
        """

rule suppa2_generateevents:
    input:
        gtf=GTF,
    output:
        ioe=pjoin(ODIR, "salmon_suppa", "suppa2.events.ioe"),
    params:
        oprefix=pjoin(ODIR, "salmon_suppa", "suppa2")
    threads: 1
    conda: "../envs/suppa2.yaml"
    log: pjoin(ODIR, "bench", "salmon_suppa", "suppa2-generateevents.time")
    shell:
        """
        /usr/bin/time -vo {log} suppa.py generateEvents -i {input.gtf} -o {params.oprefix} -e SE SS RI -f ioe
        awk 'FNR==1 && NR!=1 {{ while (/^seqname/) getline; }} 1 {{print}}' {params.oprefix}*ioe > {output.ioe}
        """

rule suppa2:
    input:
        ioe=pjoin(ODIR, "salmon_suppa", "suppa2.events.ioe"),
        quants=expand(pjoin(ODIR, "salmon_suppa", "{sample}.salmon.k{{k}}/quant.sf"),
                      sample=FQs.keys()),
    output:
        outd=directory(pjoin(ODIR, "salmon_suppa", "suppa2.k{k}/")),
    params:
        quants=pjoin(ODIR, "salmon_suppa", "*.salmon.k{k}/quant.sf")
    threads: 1
    conda: "../envs/suppa2.yaml"
    log: pjoin(ODIR, "bench", "salmon_suppa", "suppa2.k{k}.time")
    shell:
        """
        /usr/bin/time -vo {log} bash workflow/scripts/run_suppa2.sh {input.ioe} "{params.quants}" {output.outd}
        """
