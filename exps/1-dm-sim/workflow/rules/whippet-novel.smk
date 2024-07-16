rule merge_bams:
    input:
        pjoin(
            ODIR, "{sample}", "novel", "STAR", "sample_1.Aligned.sortedByCoord.out.bam"
        ),
        pjoin(
            ODIR, "{sample}", "novel", "STAR", "sample_2.Aligned.sortedByCoord.out.bam"
        ),
    output:
        pjoin(ODIR, "{sample}", "novel", "STAR", "both.bam"),
    params:
        pjoin(ODIR, "{sample}", "novel", "STAR", "both.unsrt.bam"),
    conda:
        pjoin(ENVS, "samtools.yaml")
    shell:
        """
        samtools merge {params} {input}
        samtools index {params}
        samtools sort {params} | samtools rmdup -S - {output}
        samtools index {output}
        """


rule whippet_index_novel:
    input:
        exe_i=rules.download_whippet.output.exe_i,
        fa=FA,
        gtf=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf"),
        bam=pjoin(ODIR, "{sample}", "novel", "STAR", "both.bam"),
    output:
        jls=pjoin(ODIR, "{sample}", "novel", "whippet", "index.jls"),
    params:
        index_prefix=pjoin(ODIR, "{sample}", "novel", "whippet", "index"),
    log:
        pjoin(ODIR, "{sample}", "bench", "novel", "whippet", "index.time"),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe_i} --fasta {input.fa} --gtf {input.gtf} --index {params.index_prefix} --bam {input.bam}
        """


rule whippet_quant_novel:
    input:
        exe_q=rules.download_whippet.output.exe_q,
        jls=pjoin(ODIR, "{sample}", "novel", "whippet", "index.jls"),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    output:
        pjoin(ODIR, "{sample}", "novel", "whippet", "output_{x}.psi.gz"),
    params:
        index_prefix=pjoin(ODIR, "{sample}", "novel", "whippet", "index"),
        output_prefix=pjoin(ODIR, "{sample}", "novel", "whippet", "output_{x}"),
    log:
        pjoin(ODIR, "{sample}", "bench", "novel", "whippet", "quant-{x}.time"),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe_q} --index {params.index_prefix} --out {params.output_prefix} --biascorrect {input.fq1} {input.fq2}
        """


rule whippet_delta_novel:
    input:
        exe_d=rules.download_whippet.output.exe_d,
        psi1=pjoin(ODIR, "{sample}", "novel", "whippet", "output_1.psi.gz"),
        psi2=pjoin(ODIR, "{sample}", "novel", "whippet", "output_2.psi.gz"),
    output:
        gz=pjoin(ODIR, "{sample}", "novel", "whippet", "psi.diff.gz"),
        diff=pjoin(ODIR, "{sample}", "novel", "whippet", "psi.diff"),
    params:
        prefix=pjoin(ODIR, "{sample}", "novel", "whippet", "psi"),
    log:
        pjoin(ODIR, "{sample}", "bench", "novel", "whippet", "delta.time"),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe_d} -a {input.psi1}, -b {input.psi2}, -o {params.prefix}
        gunzip -k {output.gz}
        """
