rule rmats_c1:
    input:
        pjoin(
            ODIR,
            "{sample}",
            "{annov}",
            "STAR",
            "sample_1.Aligned.sortedByCoord.out.bam",
        ),
    output:
        pjoin(ODIR, "{sample}", "{annov}", "STAR.1.txt"),
    shell:
        """
        echo {input} > {output}
        """


rule rmats_c2:
    input:
        pjoin(
            ODIR,
            "{sample}",
            "{annov}",
            "STAR",
            "sample_2.Aligned.sortedByCoord.out.bam",
        ),
    output:
        pjoin(ODIR, "{sample}", "{annov}", "STAR.2.txt"),
    shell:
        """
        echo {input} > {output}
        """


rule rmats:
    input:
        gtf=lambda wildcards: (
            pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf")
            if wildcards.annov == "anno"
            else pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf")
        ),
        c1txt=pjoin(ODIR, "{sample}", "{annov}", "STAR.1.txt"),
        c2txt=pjoin(ODIR, "{sample}", "{annov}", "STAR.2.txt"),
    output:
        outd=directory(pjoin(ODIR, "{sample}", "{annov}", "rMATS")),
        summary=pjoin(ODIR, "{sample}", "{annov}", "rMATS", "summary.txt"),
    params:
        tmpd=pjoin(ODIR, "{sample}", "{annov}", "rMATS-tmp"),
        novelp=lambda wildcards: (
            "--novelSS --mil 30 --mel 5000" if wildcards.annov == "novel" else ""
        ),
    threads: workflow.cores
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "rMATS.time"),
    conda:
        pjoin(ENVS, "rmats.yaml")
    shell:
        """
        /usr/bin/time -vo {log} rmats.py {params.novelp} --gtf {input.gtf} --b1 {input.c1txt} --b2 {input.c2txt} --od {output.outd} --tmp {params.tmpd} --readLength 150 --nthread {threads} -t paired
        """


rule build_rmats:
    input:
        truth=pjoin(ODIR, "{sample}", "truth.csv"),
        calls=pjoin(ODIR, "{sample}", "{annov}", "rMATS"),
    output:
        pjoin(ODIR, "{sample}", "rMATS-{annov}.csv"),
    shell:
        """
        python3 build_rmats.py {wildcards.annov} {input.calls} > {output}
        """
