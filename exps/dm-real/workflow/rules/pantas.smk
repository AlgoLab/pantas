rule download_pantas:
    output:
        exe=pjoin(software_folder, "pantas2", "pantas2.sh"),
        augm=pjoin(software_folder, "pantas2", "scripts", "alignments_augmentation_from_gaf.py"),
        call=pjoin(software_folder, "pantas2", "scripts", "call.py"),
        quant=pjoin(software_folder, "pantas2", "scripts", "quantify3.py"),
        outd=directory(pjoin(software_folder, "pantas2")),
    shell:
        """
        rm -r {output.outd}
        git clone git@github.com:AlgoLab/pantas2.git {output.outd}
        """


rule pantas2_index:
    input:
        fa=FA,
        gtf=GTF,
        vcf=VCF,
        exe=pjoin(software_folder, "pantas2", "pantas2.sh"),
    output:
        xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.dist"),
        gfa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.annotated.gfa"),
        splicedfa=pjoin(ODIR, "pantas2", "index", "genes.spliced.fa"),
    params:
        wd=pjoin(ODIR,  "pantas2", "index"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "index.time"),
    threads: workflow.cores
    conda: "../envs/pantas2.yaml"
    shell:
        """
        curr_dir=$(pwd)
        /usr/bin/time -vo {log.time} bash {input.exe} index {input.fa} {input.gtf} {input.vcf} {params.wd} {threads}
        """



rule pantas2_mpmap:
    input:
        xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.dist"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
    output:
        gaf=pjoin(ODIR,  "pantas2", "{sample}.gaf"),
    threads: workflow.cores
    log:
        pjoin(ODIR, "bench", "pantas2", "mpmap{sample}.time"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log} vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F GAF --threads {threads} > {output.gaf}
        """

rule pantas_weight:
    input:
        gfa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.annotated.gfa"),
        gaf=pjoin(ODIR,  "pantas2", "{sample}.gaf"),
        augm=pjoin(software_folder, "pantas2", "scripts", "alignments_augmentation_from_gaf.py"),
    output:
        gfa=pjoin(ODIR,  "pantas2", "graph_{sample}.gfa"),
    log:
        pjoin(ODIR, "bench", "pantas2", "weight{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 {input.augm} {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call:
    input:
        gfa=pjoin(ODIR,  "pantas2", "graph_{sample}.gfa"),
        gtf=GTF,
        call=pjoin(software_folder, "pantas2", "scripts", "call.py"),
    output:
        csv=pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"),
    log:
        pjoin(ODIR,  "bench",  "pantas2", "call{sample}.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 {input.call} --rc {wildcards.w} {input.gfa} {input.gtf} > {output.csv}
        """



rule c1_csv:
    input:
        expand(pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"), sample
               = C1.keys(), w="{w}"),
    output:
        pjoin(ODIR, "pantas2", "{w}", "c1.csv"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        csvstack {input} > {output}
        """


rule c2_csv:
    input:
        expand(pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"), sample
               = C2.keys(), w="{w}"),
    output:
        pjoin(ODIR, "pantas2", "{w}", "c2.csv"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        csvstack {input} > {output}
        """

rule pantas_quant:
    input:
        csv1s=expand(pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"), sample
               = C1.keys(), w="{w}"),
        csv2s=expand(pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"), sample
               = C2.keys(), w="{w}"),
        quant=pjoin(software_folder, "pantas2", "scripts", "quantify3.py"),
    output:
        csv=pjoin(ODIR, "pantas2", "quant.w{w}.csv"),
    log:
        pjoin(ODIR, "bench", "pantas2", "quant.w{w}.time"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log} python3 {input.quant} -c1 {input.csv1s} -c2 {input.csv2s} > {output.csv}
        """
