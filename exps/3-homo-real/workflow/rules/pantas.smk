rule sub_fa:
    input:
        fa=FA,
        bed=pjoin(ODIR, "pantas2", "genes.bed"),
    output:
        fa=pjoin(ODIR, "pantas2", "ref.fa"),
    shell:
        """
        cut -f1 {input.bed} | sort -u | while read chr ; do samtools faidx {input.fa} ${{chr}} ; done > {output.fa}
        samtools faidx {output.fa}
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
    threads: 16
    log:
        time=pjoin(ODIR, "bench", "pantas2", "shark_{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} shark --threads {threads} -q 10 -r {input.fa} -1 {input.fq1} -2 {input.fq2} -o {output.fq1} -p {output.fq2} > {output.tsv}
        """


rule pantas2_index:
    input:
        fa=pjoin(ODIR, "pantas2", "ref.fa"),
        gtf=pjoin(ODIR, "pantas2", "genes.gtf"),
        vcf=pjoin(ODIR, "pantas2", "variations.vcf.gz"),
    output:
        xg=pjoin(ODIR, "pantas2", "index", "pantranscriptome.xg"),
        gfa=pjoin(ODIR, "pantas2", "index", "pantranscriptome-annotated.gfa"),
    params:
        wd=pjoin(ODIR, "pantas2", "index"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "build.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../../pantas build -r -t {threads} -o {params.wd} {input.fa} {input.gtf} {input.vcf}
        """

rule gcsa2_index:
    input:
        xg=pjoin(ODIR, "pantas2", "index", "pantranscriptome.xg"),
    output:
        gcsa=pjoin(
            ODIR, "pantas2", "index", "pantranscriptome.gcsa"
        ),
        dist=pjoin(
            ODIR, "pantas2", "index", "pantranscriptome.dist"
        ),
    params:
        wd=pjoin(ODIR, "pantas2", "index"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "index.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} vg index --progress --threads {threads} --temp-dir {params.wd} --gcsa-out {output.gcsa} --dist-name {output.dist} {input.xg}
        """

# rule get_intron_distr:
#     input:
#         gtf=pjoin(ODIR, "pantas2", "genes.gtf"),
#     output:
#         distr=pjoin(ODIR, "pantas2", "introns.distr"),
#     conda:
#         "../envs/pylib.yaml"
#     shell:
#         """
#         python3 workflow/scripts/intron_length_distribution.py -g {input.gtf} -o {output.distr}
#         """


rule pantas2_mpmap:
    input:
        xg=pjoin(ODIR, "pantas2", "index", "pantranscriptome.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "pantranscriptome.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "pantranscriptome.dist"),
        fq1=pjoin(ODIR, "pantas2", "{sample}_1.fq"),
        fq2=pjoin(ODIR, "pantas2", "{sample}_2.fq"),
        # distr=pjoin(ODIR, "pantas2", "introns.distr"),
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
# --intron-distr {input.distr}

rule pantas_weight:
    input:
        gfa=pjoin(ODIR, "pantas2", "index", "pantranscriptome-annotated.gfa"),
        gaf=pjoin(ODIR, "pantas2", "{sample}.gaf"),
    output:
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "weight-{sample}.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../../pantas augment {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call:
    input:
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
        gtf=rules.sub_gtf.output.gtf,
    output:
        csv=pjoin(ODIR, "pantas2", "w{w}", "events_{sample}.csv"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "call-{sample}.w{w}.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../../pantas call -e ES -n -w {wildcards.w} {input.gfa} {input.gtf} > {output.csv}
        """


rule pantas_quant:
    input:
        csv1=expand(
            pjoin(ODIR, "pantas2", "w{w}", "events_{sample}.csv"),
            sample=C1.keys(),
            w="{w}",
        ),
        csv2=expand(
            pjoin(ODIR, "pantas2", "w{w}", "events_{sample}.csv"),
            sample=C2.keys(),
            w="{w}",
        ),
    output:
        csv=pjoin(ODIR, "pantas2", "quant.w{w}.csv"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "quant.w{w}.time"),
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../../pantas quant -a {input.csv1} {input.csv2} > {output.csv}
        """

rule pantas_remap:
    input:
        csv=pjoin(ODIR, "pantas2", "quant.w{w}.csv"),
        gtf=rules.sub_gtf.output.gtf,
    output:
        csv=pjoin(ODIR, "pantas2", "quant-remap.w{w}.csv"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "remap.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../../pantas remap {input.csv} {input.gtf} > {output.csv}
        """
