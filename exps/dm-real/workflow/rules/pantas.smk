rule download_pantas:
    output:
        exe=pjoin(software_folder, "pantas", "index.smk"),
        augm=pjoin(software_folder, "pantas", "scripts", "alignments_augmentation_from_gaf.py"),
        call=pjoin(software_folder, "pantas", "scripts", "call.py"),
        quant=pjoin(software_folder, "pantas", "scripts", "quantify3.py"),
        outd=directory(pjoin(software_folder, "pantas")),
    shell:
        """
        rm -r {output.outd}
        git clone git@github.com:AlgoLab/pantas.git {output.outd}
        """


rule pantas2_index:
    input:
        fa=FA,
        gtf=GTF,
        vcf=VCF,
        exe=pjoin(software_folder, "pantas", "index.smk"),
    output:
        xg=pjoin(ODIR, "pantas", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, "pantas", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, "pantas", "index", "spliced-pangenome.dist"),
        gfa=pjoin(ODIR, "pantas", "index", "spliced-pangenome.annotated.gfa"),
        splicedfa=pjoin(ODIR, "pantas", "index", "genes.spliced.fa"),
    params:
        wd=pjoin(ODIR, "pantas", "index"),
    log:
        time=pjoin(ODIR, "bench", "pantas", "index.time"),
    threads: workflow.cores
    conda: "../envs/pantas2.yaml"
    shell:
        """
        pushd {software_folder}/pantas
        /usr/bin/time -vo {log.time} snakemake -s index.smk -c{threads} --config fa={input.fa} gtf={input.gtf} vcf={input.vcf} wd={params.wd}
        popd
        """



rule pantas2_mpmap:
    input:
        xg=pjoin(ODIR, "pantas", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, "pantas", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, "pantas", "index", "spliced-pangenome.dist"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
    output:
        gaf=pjoin(ODIR,  "pantas", "{sample}.gaf"),
    threads: workflow.cores
    log:
        pjoin(ODIR, "bench", "pantas", "mpmap{sample}.time"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log} vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F GAF --threads {threads} > {output.gaf}
        """

rule pantas_weight:
    input:
        gfa=pjoin(ODIR, "pantas", "index", "spliced-pangenome.annotated.gfa"),
        gaf=pjoin(ODIR,  "pantas", "{sample}.gaf"),
        augm=pjoin(software_folder, "pantas", "scripts", "alignments_augmentation_from_gaf.py"),
    output:
        gfa=pjoin(ODIR,  "pantas", "graph_{sample}.gfa"),
    log:
        pjoin(ODIR, "bench", "pantas", "weight{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 {input.augm} {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call:
    input:
        gfa=pjoin(ODIR,  "pantas", "graph_{sample}.gfa"),
        gtf=GTF,
        call=pjoin(software_folder, "pantas", "scripts", "call.py"),
    output:
        csv=pjoin(ODIR, "pantas", "{w}", "events_{sample}.csv"),
    log:
        pjoin(ODIR,  "bench",  "pantas", "call{sample}.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 {input.call} --rca {wildcards.w} --rc {wildcards.w} {input.gfa} {input.gtf} > {output.csv}
        """



rule c1_csv:
    input:
        expand(pjoin(ODIR, "pantas", "{w}", "events_{sample}.csv"), sample
               = C1.keys(), w="{w}"),
    output:
        pjoin(ODIR, "pantas", "{w}", "c1.csv"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        csvstack {input} > {output}
        """


rule c2_csv:
    input:
        expand(pjoin(ODIR, "pantas", "{w}", "events_{sample}.csv"), sample
               = C2.keys(), w="{w}"),
    output:
        pjoin(ODIR, "pantas", "{w}", "c2.csv"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        csvstack {input} > {output}
        """

rule pantas_quant:
    input:
        csv1s=expand(pjoin(ODIR, "pantas", "{w}", "events_{sample}.csv"), sample
               = C1.keys(), w="{w}"),
        csv2s=expand(pjoin(ODIR, "pantas", "{w}", "events_{sample}.csv"), sample
               = C2.keys(), w="{w}"),
        quant=pjoin(software_folder, "pantas", "scripts", "quantify3.py"),
    output:
        csv=pjoin(ODIR, "pantas", "quant.w{w}.csv"),
    log:
        pjoin(ODIR, "bench", "pantas", "quant.w{w}.time"),
    conda: "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log} python3 {input.quant} -c1 {input.csv1s} -c2 {input.csv2s} > {output.csv}
        """
