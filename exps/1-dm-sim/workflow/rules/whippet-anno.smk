rule download_whippet:
    output:
        exe_i=pjoin(ODIR, "software", "whippet", "bin", "whippet-index.jl"),
        exe_q=pjoin(ODIR, "software", "whippet", "bin", "whippet-quant.jl"),
        exe_d=pjoin(ODIR, "software", "whippet", "bin", "whippet-delta.jl"),
        outd=directory(pjoin(ODIR, "software", "whippet")),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        rm -r {output.outd}
        git clone https://github.com/timbitz/Whippet.jl.git {output.outd}
        cd {output.outd}
        julia --project -e \'using Pkg; Pkg.instantiate(); \'
        """


rule whippet_index_anno:
    input:
        exe_i=rules.download_whippet.output.exe_i,
        fa=FA,
        gtf=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf"),
    output:
        jls=pjoin(ODIR, "{sample}", "anno", "whippet", "index.jls"),
    params:
        index_prefix=pjoin(ODIR, "{sample}", "anno", "whippet", "index"),
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "whippet", "index.time"),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe_i} --fasta {input.fa} --gtf {input.gtf} --index {params.index_prefix}
        """


rule whippet_quant_anno:
    input:
        exe_q=rules.download_whippet.output.exe_q,
        jls=pjoin(ODIR, "{sample}", "anno", "whippet", "index.jls"),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    output:
        pjoin(ODIR, "{sample}", "anno", "whippet", "output_{x}.psi.gz"),
    params:
        index_prefix=pjoin(ODIR, "{sample}", "anno", "whippet", "index"),
        output_prefix=pjoin(ODIR, "{sample}", "anno", "whippet", "output_{x}"),
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "whippet", "quant-{x}.time"),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe_q} --index {params.index_prefix} --out {params.output_prefix} --biascorrect {input.fq1} {input.fq2} 
        """


rule whippet_delta_anno:
    input:
        exe_d=rules.download_whippet.output.exe_d,
        psi1=pjoin(ODIR, "{sample}", "anno", "whippet", "output_1.psi.gz"),
        psi2=pjoin(ODIR, "{sample}", "anno", "whippet", "output_2.psi.gz"),
    output:
        gz=pjoin(ODIR, "{sample}", "anno", "whippet", "psi.diff.gz"),
        diff=pjoin(ODIR, "{sample}", "anno", "whippet", "psi.diff"),
    params:
        prefix=pjoin(ODIR, "{sample}", "anno", "whippet", "psi"),
    log:
        pjoin(ODIR, "{sample}", "bench", "anno", "whippet", "delta.time"),
    conda:
        pjoin(ENVS, "julia.yaml")
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe_d} -a {input.psi1}, -b {input.psi2}, -o {params.prefix}
        gunzip -k {output.gz}
        """
