rule download_whippet:
    output:
        exe_i=pjoin(software_folder, "whippet", "bin", "whippet-index.jl"),
        exe_q=pjoin(software_folder, "whippet", "bin", "whippet-quant.jl"),
        exe_d=pjoin(software_folder, "whippet", "bin", "whippet-delta.jl"),
        outd=directory(pjoin(software_folder, "whippet")),
    conda:
        "../envs/julia.yaml"
    shell:
        """
        rm -r {output.outd}
        git clone https://github.com/timbitz/Whippet.jl.git {output.outd}
        cd {output.outd}
        julia --project -e \'using Pkg; Pkg.instantiate(); \'
        """


rule whippet_index:
    input:
        fa=FA,
        gtf=GTF,
        exe=pjoin(software_folder, "whippet", "bin", "whippet-index.jl"),
    output:
        jls=pjoin(ODIR, "whippet", "index.jls"),
    params:
        index_prefix=pjoin(ODIR, "whippet", "index"),
    conda:
        "../envs/julia.yaml"
    log:
        pjoin(ODIR, "bench", "whippet", "index.time"),
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe} --fasta {input.fa} --gtf {input.gtf} --index {params.index_prefix} --suppress-low-tsl
        """


rule whippet_quant:
    input:
        jls=pjoin(ODIR, "whippet", "index.jls"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
        exe=pjoin(software_folder, "whippet", "bin", "whippet-quant.jl"),
    output:
        pjoin(ODIR, "whippet", "{sample}.psi.gz"),
    params:
        index_prefix=pjoin(ODIR, "whippet", "index"),
        output_prefix=pjoin(ODIR, "whippet", "{sample}"),
    conda:
        "../envs/julia.yaml"
    log:
        pjoin(ODIR, "bench", "whippet", "quant-{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log} julia {input.exe} --index {params.index_prefix} --out {params.output_prefix} --biascorrect {input.fq1} {input.fq2}
        """


rule whippet_c1:
    input:
        expand(pjoin(ODIR, "whippet", "{sample}.psi.gz"), sample=C1.keys()),
    output:
        pjoin(ODIR, "whippet", "c1-psis.txt"),
    shell:
        """
        echo {input} | tr " " "," > {output}
        """


rule whippet_c2:
    input:
        expand(pjoin(ODIR, "whippet", "{sample}.psi.gz"), sample=C2.keys()),
    output:
        pjoin(ODIR, "whippet", "c2-psis.txt"),
    shell:
        """
        echo {input} | tr " " "," > {output}
        """


rule whippet_delta:
    input:
        psi1=pjoin(ODIR, "whippet", "c1-psis.txt"),
        psi2=pjoin(ODIR, "whippet", "c2-psis.txt"),
        exe=pjoin(software_folder, "whippet", "bin", "whippet-delta.jl"),
    output:
        gz=pjoin(ODIR, "whippet", "psi.diff.gz"),
        diff=pjoin(ODIR, "whippet", "psi.diff"),
    params:
        prefix=pjoin(ODIR, "whippet", "psi"),
    conda:
        "../envs/julia.yaml"
    log:
        pjoin(ODIR, "bench", "whippet", "delta.time"),
    shell:
        """
        A=$(cat {input.psi1})
        B=$(cat {input.psi2})
        /usr/bin/time -vo {log} julia {input.exe} -a $A -b $B -o {params.prefix}
        gunzip -k {output.gz}
        """
