rule suppa2_generateevents:
    input:
        gtf=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf"),
    output:
        ioe=pjoin(ODIR, "{sample}", "anno", "suppa2", "annotated.events.ioe"),
    params:
        oprefix=pjoin(ODIR, "{sample}", "anno", "suppa2", "annotated"),
    threads: 1
    conda:
        pjoin(ENVS, "suppa2.yaml")
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "suppa2", "generateevents.time"),
    shell:
        """
        /usr/bin/time -vo {log} suppa.py generateEvents -i {input.gtf} -o {params.oprefix} -e SE SS RI -f ioe
        awk 'FNR==1 && NR!=1 {{ while (/^seqname/) getline; }} 1 {{print}}' {params.oprefix}*ioe > {output.ioe}
        """


rule suppa2:
    input:
        ioe=pjoin(ODIR, "{sample}", "anno", "suppa2", "annotated.events.ioe"),
        quant1=pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-quant-1", "quant.sf"),
        quant2=pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-quant-2", "quant.sf"),
    output:
        outd=directory(pjoin(ODIR, "{sample}", "anno", "suppa2", "OUT")),
        psi=pjoin(ODIR, "{sample}", "anno", "suppa2", "OUT", "DIFF.dpsi"),
    params:
        quants=pjoin(ODIR, "{sample}", "anno", "salmon", "salmon-quant-*", "quant.sf"),
    threads: 1
    conda:
        pjoin(ENVS, "suppa2.yaml")
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "suppa2", "suppa2.time"),
    shell:
        """
        # NOTE: script is update to work with single replicate per condition
        /usr/bin/time -vo {log} bash run_suppa2.sh {input.ioe} "{params.quants}" {output.outd}
        """


rule build_suppa:
    input:
        calls=pjoin(ODIR, "{sample}", "anno", "suppa2", "OUT", "DIFF.dpsi"),
    output:
        csv=pjoin(ODIR, "{sample}", "SUPPA2-anno.csv"),
    shell:
        """
        python3 build_suppa.py {input.calls} 100 > {output.csv}
        """
