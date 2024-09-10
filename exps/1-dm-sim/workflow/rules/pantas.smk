rule remove_sample_from_vcf:
    input:
        vcf=VCF,
    output:
        vcf=pjoin(ODIR, "{sample}", "population-oneout.vcf.gz"),
    shell:
        """
        bcftools view -s ^{wildcards.sample} -Oz {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule pantas2_build:
    input:
        fa=FA,
        gtf=lambda wildcards: (
            pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf")
            if wildcards.annov == "anno"
            else pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf")
        ),
        vcf=pjoin(ODIR, "{sample}", "population-oneout.vcf.gz"),
    output:
        xg=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.xg"),
        gfa=pjoin(
            ODIR,
            "{sample}",
            "{annov}",
            "pantas2",
            "index",
            "pantranscriptome-annotated.gfa",
        ),
    params:
        wd=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "index"),
    log:
        time=pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "build.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../../pantas build -t {threads} -o {params.wd} {input.fa} {input.gtf} {input.vcf}
        """


rule gcsa2_index:
    input:
        xg=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.xg"),
    output:
        gcsa=pjoin(
            ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.gcsa"
        ),
        dist=pjoin(
            ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.dist"
        ),
    params:
        wd=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "index"),
    log:
        time=pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "index.time"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} vg index --progress --threads {threads} --temp-dir {params.wd} --gcsa-out {output.gcsa} --dist-name {output.dist} {input.xg}
        """


rule pantas2_mpmap:
    input:
        xg=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.xg"),
        gcsa=pjoin(
            ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.gcsa"
        ),
        dist=pjoin(
            ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.dist"
        ),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    output:
        gaf=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "sample_{x}.gaf"),
    threads: workflow.cores
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "mpmap{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F GAF --threads {threads} > {output.gaf}
        """


rule pantas2_mpmap_bam:
    input:
        xg=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.xg"),
        gcsa=pjoin(
            ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.gcsa"
        ),
        dist=pjoin(
            ODIR, "{sample}", "{annov}", "pantas2", "index", "pantranscriptome.dist"
        ),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    output:
        bam=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "sample_{x}.bam"),
    threads: workflow.cores
    shell:
        """
        vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F BAM --threads {threads} | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule pantas_weight:
    input:
        gfa=pjoin(
            ODIR,
            "{sample}",
            "{annov}",
            "pantas2",
            "index",
            "pantranscriptome-annotated.gfa",
        ),
        gaf=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "sample_{x}.gaf"),
    output:
        gfa=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "graph_{x}.gfa"),
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "weight{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} bash ../../pantas augment {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call:
    input:
        gfa=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "graph_{x}.gfa"),
        gtf=lambda wildcards: (
            pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf")
            if wildcards.annov == "anno"
            else pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf")
        ),
    output:
        csv=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "events_{x}.w{w}.csv"),
        # log=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "events_{x}.w{w}.log"),
    params:
        novelp=lambda wildcards: "" if wildcards.annov == "anno" else "-n -a",
    threads: workflow.cores / 2
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "call{x}.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log} bash ../../pantas call -l 25 {params.novelp} -w {wildcards.w} {input.gfa} {input.gtf} > {output.csv}
        """


rule pantas_quant:
    input:
        csv1=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "events_1.w{w}.csv"),
        csv2=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "events_2.w{w}.csv"),
    output:
        csv=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "quant.w{w}.csv"),
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "quant.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log} bash ../../pantas quant {input.csv1} {input.csv2} > {output.csv}
        """


rule pantas_remap:
    input:
        csv=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "quant.w{w}.csv"),
        gtf=lambda wildcards: (
            pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf")
            if wildcards.annov == "anno"
            else pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf")
        ),
    output:
        csv=pjoin(ODIR, "{sample}", "{annov}", "pantas2", "quant-remap.w{w}.csv"),
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "pantas2", "remap.w{w}.time"),
    shell:
        """
        /usr/bin/time -vo {log} bash ../../pantas remap -i 25 {input.csv} {input.gtf} > {output.csv}
        """
