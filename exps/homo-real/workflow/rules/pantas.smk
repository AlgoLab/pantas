rule download_pantas:
    output:
        sh=pjoin(software_folder, "pantas2", "pantas2.sh"),
        augm=pjoin(
            software_folder,
            "pantas2",
            "scripts",
            "alignments_augmentation_from_gaf.py",
        ),
        call=pjoin(software_folder, "pantas2", "scripts", "call.py"),
        quant=pjoin(software_folder, "pantas2", "scripts", "quantify.py"),
        outd=directory(pjoin(software_folder, "pantas2")),
    shell:
        """
        rm -r {output.outd}
        git clone git@github.com:AlgoLab/pantas2.git {output.outd}
        """


rule sub_fa:
    input:
        fa=FA,
        bed=pjoin(ODIR, "pantas2", "genes.bed"),
    output:
        fa=pjoin(ODIR, "pantas2", "ref.fa"),
    shell:
        """
        cut -f1 {input.bed} | sort -u | while read chr ; do samtools faidx {input.fa} ${{chr}} ; done > {output.fa}
        """


rule sub_gtf:
    input:
        gtf=GTF,
        genes=GENES,
    output:
        gtf=pjoin(ODIR, "pantas2", "genes.gtf"),
        bed=pjoin(ODIR, "pantas2", "genes.bed"),
    shell:
        """
        grep -f {input.genes} {input.gtf} > {output.gtf}
        grep -P "\\tgene\\t" {output.gtf} | cut -f1,4,5 > {output.bed}
        """


rule sub_vcf:
    input:
        vcf=VCF,
        bed=rules.sub_gtf.output.bed,
    output:
        vcf=pjoin(ODIR, "pantas2", "variations.vcf.gz"),
    conda:
        "../envs/biotools.yaml"
    shell:
        """        
        bcftools view -Oz -R {input.bed} {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule get_genes_fa:
    input:
        fa=FA,
        gtf=rules.sub_gtf.output.gtf,
    output:
        fa=pjoin(ODIR, "pantas2", "genes.fa"),
    conda:
        "../envs/biotools.yaml"
    shell:
        """
        bash workflow/scripts/get_genes_fa.sh {input.fa} {input.gtf} > {output.fa}
        """


rule shark:
    input:
        fa=pjoin(ODIR, "pantas2", "genes.fa"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
    output:
        fq1=pjoin(ODIR, "pantas2", "{sample}_1.fq"),
        fq2=pjoin(ODIR, "pantas2", "{sample}_2.fq"),
        tsv=pjoin(ODIR, "pantas2", "{sample}.tsv"),
    conda:
        "../envs/shark.yaml"
    log:
        time=pjoin(ODIR, "bench", "pantas2", "shark_{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} shark -r {input.fa} -1 {input.fq1} -2 {input.fq2} -o {output.fq1} -p {output.fq2} > {output.tsv}
        """


rule pantas2_index:
    input:
        fa=pjoin(ODIR, "pantas2", "ref.fa"),
        gtf=pjoin(ODIR, "pantas2", "genes.gtf"),
        vcf=pjoin(ODIR, "pantas2", "variations.vcf.gz"),
        exe=rules.download_pantas.output.sh,
    output:
        xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.dist"),
        gfa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.annotated.gfa"),
        splicedfa=pjoin(ODIR, "pantas2", "index", "genes.spliced.fa"),
    params:
        wd=pjoin(ODIR, "pantas2", "index"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "index.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} bash {input.exe} index {input.fa} {input.gtf} {input.vcf} {params.wd} {threads}
        """


rule pantas2_mpmap:
    input:
        xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.dist"),
        fq1=pjoin(ODIR, "pantas2", "{sample}_1.fq"),
        fq2=pjoin(ODIR, "pantas2", "{sample}_2.fq"),
    output:
        gaf=pjoin(ODIR, "pantas2", "{sample}.gaf"),
    threads: workflow.cores
    log:
        time=pjoin(ODIR, "bench", "pantas2", "mpmap-{sample}.time"),
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F GAF --threads {threads} > {output.gaf}
        """


# rule pantas2_mpmap_bam:
#     input:
#         xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.xg"),
#         gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.gcsa"),
#         dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.dist"),
#         fq1=pjoin(ODIR, "pantas2", "{sample}_1.fq"),
#         fq2=pjoin(ODIR, "pantas2", "{sample}_2.fq"),
#     output:
#         bam=pjoin(ODIR, "pantas2", "{sample}.bam"),
#     threads: workflow.cores
#     conda:
#         "../envs/pantas2.yaml"
#     shell:
#         """
#         vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F BAM --threads {threads} | samtools sort > {output.bam}
#         samtools index {output.bam}
#         """


rule pantas_weight:
    input:
        gfa=pjoin(ODIR, "pantas2", "index", "spliced-pangenome.annotated.gfa"),
        gaf=pjoin(ODIR, "pantas2", "{sample}.gaf"),
        exe=rules.download_pantas.output.augm,
    output:
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "weight-{sample}.time"),
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} python3 {input.exe} {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call:
    input:
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
        gtf=rules.sub_gtf.output.gtf,
        exe=rules.download_pantas.output.call,
    output:
        csv=pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "call-{sample}.w{w}.time"),
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} python3 {input.exe} --rc {wildcards.w} {input.gfa} {input.gtf} > {output.csv}
        """


# rule c1_csv:
#     input:
#         expand(
#             pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"),
#             sample=C1.keys(),
#             w="{w}",
#         ),
#     output:
#         pjoin(ODIR, "pantas2", "{w}", "c1.csv"),
#     conda:
#         "../envs/pantas2.yaml"
#     shell:
#         """
#         csvstack {input} > {output}
#         """
# rule c2_csv:
#     input:
#         expand(
#             pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"),
#             sample=C2.keys(),
#             w="{w}",
#         ),
#     output:
#         pjoin(ODIR, "pantas2", "{w}", "c2.csv"),
#     conda:
#         "../envs/pantas2.yaml"
#     shell:
#         """
#         csvstack {input} > {output}
#         """


rule pantas_quant:
    input:
        splicedfa=pjoin(ODIR, "pantas2", "index", "genes.spliced.fa"),
        csv1s=expand(
            pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"),
            sample=C1.keys(),
            w="{w}",
        ),
        csv2s=expand(
            pjoin(ODIR, "pantas2", "{w}", "events_{sample}.csv"),
            sample=C2.keys(),
            w="{w}",
        ),
        exe=rules.download_pantas.output.quant,
    output:
        csv=pjoin(ODIR, "pantas2", "quant.w{w}.csv"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "quant.w{w}.time"),
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} python3 {input.exe} {input.splicedfa} {input.csv1s} {input.csv2s} > {output.csv}
        """
