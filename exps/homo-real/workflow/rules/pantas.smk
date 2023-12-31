rule download_pantas:
    output:
        index=pjoin(software_folder, "pantas2", "index-reduced.smk"),
        augm=pjoin(
            software_folder,
            "pantas2",
            "scripts",
            "alignments_augmentation_from_gaf.py",
        ),
        call=pjoin(software_folder, "pantas2", "scripts", "call.py"),
        quant=pjoin(software_folder, "pantas2", "scripts", "quantify3.py"),
        outd=directory(pjoin(software_folder, "pantas2")),
    shell:
        """
        rm -r {output.outd}
        git clone git@github.com:AlgoLab/pantas.git {output.outd}
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
        smk=rules.download_pantas.output.index,
    output:
        xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.dist"),
        gfa=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.annotated.gfa"),
        splicedfa=pjoin(ODIR, "pantas2", "index", "genes.spliced.fa"),
        gfarp=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.refpath"),
    params:
        wd=pjoin(ODIR, "pantas2", "index"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "index.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: workflow.cores
    shell:
        """
        pushd {software_folder}/pantas2
        /usr/bin/time -vo {log.time} snakemake -s index-reduced.smk -c{threads} --config fa={input.fa} gtf={input.gtf} vcf={input.vcf} wd={params.wd}
        popd
        """


rule get_intron_distr:
    input:
        gtf=pjoin(ODIR, "pantas2", "genes.gtf"),
    output:
        distr=pjoin(ODIR, "pantas2", "introns.distr"),
    conda:
        "../envs/pylib.yaml"
    shell:
        """
        python3 workflow/scripts/intron_length_distribution.py -g {input.gtf} -o {output.distr}
        """


rule pantas2_mpmap:
    input:
        xg=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.xg"),
        gcsa=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.gcsa"),
        dist=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.dist"),
        fq1=pjoin(ODIR, "pantas2", "{sample}_1.fq"),
        fq2=pjoin(ODIR, "pantas2", "{sample}_2.fq"),
        distr=pjoin(ODIR, "pantas2", "introns.distr"),
    output:
        gaf=pjoin(ODIR, "pantas2", "{sample}.gaf"),
    threads: 16  # workflow.cores
    log:
        time=pjoin(ODIR, "bench", "pantas2", "mpmap-{sample}.time"),
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} vg mpmap --intron-distr {input.distr} -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F GAF --threads {threads} > {output.gaf}
        """


rule pantas_weight:
    input:
        gfa=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.annotated.gfa"),
        gaf=pjoin(ODIR, "pantas2", "{sample}.gaf"),
        exe=rules.download_pantas.output.augm,
    output:
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "weight-{sample}.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} python3 {input.exe} {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call:
    input:
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
        gtf=rules.sub_gtf.output.gtf,
        gfarp=rules.pantas2_index.output.gfarp,
        exe=rules.download_pantas.output.call,
    output:
        csv=pjoin(ODIR, "pantas2", "w{w}", "events_{sample}.csv"),
    log:
        time=pjoin(ODIR, "bench", "pantas2", "call-{sample}.w{w}.time"),
    conda:
        "../envs/pantas2.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} python3 {input.exe} --events ES --novel --rca -1 --rc {wildcards.w} {input.gfa} {input.gtf} --rp {input.gfarp} > {output.csv}
        """


rule pantas_quant:
    input:
        csv1s=expand(
            pjoin(ODIR, "pantas2", "w{w}", "events_{sample}.csv"),
            sample=C1.keys(),
            w="{w}",
        ),
        csv2s=expand(
            pjoin(ODIR, "pantas2", "w{w}", "events_{sample}.csv"),
            sample=C2.keys(),
            w="{w}",
        ),
        gfa1s=expand(
            pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
            sample=C1.keys(),
            w="{w}",
        ),
        gfa2s=expand(
            pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
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
        /usr/bin/time -vo {log.time} python3 {input.exe} -c1 {input.csv1s} -c2 {input.csv2s} --relax 0 --minj 10 --gfac1 {input.gfa1s} --gfac2 {input.gfa2s} > {output.csv}
        """


rule pantas2_gaf2sam:
    input:
        gaf=pjoin(ODIR, "pantas2", "{sample}.gaf"),
        gfa=pjoin(ODIR, "pantas2", "graph_{sample}.gfa"),
        gfarp=pjoin(ODIR, "pantas2", "index", "spliced-pangenes.refpath"),
    output:
        sam=pjoin(ODIR, "pantas2", "{sample}.sam"),
    log:
        err=pjoin(ODIR, "pantas2", "{sample}.err"),
    threads: 1
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        python3 workflow/scripts/gaf2sam.py {input.gaf} {input.gfa} {input.gfarp} > {output.sam} 2> {log.err}
        """


rule sam2bam:
    input:
        sam=pjoin(ODIR, "pantas2", "{sample}.sam"),
    output:
        bam=pjoin(ODIR, "pantas2", "{sample}.bam"),
    threads: 1
    conda:
        "../envs/pantas2.yaml"
    shell:
        """
        samtools view -bS {input.sam} | samtools sort > {output.bam}
        samtools index {output.bam}
        """
