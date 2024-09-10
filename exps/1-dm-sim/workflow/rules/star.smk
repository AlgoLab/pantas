rule STAR_index_anno:
    input:
        fa=FA,
        gtf=lambda wildcards: (
            pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf")
            if wildcards.annov == "anno"
            else pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf")
        ),
    output:
        index=directory(pjoin(ODIR, "{sample}", "{annov}", "STAR-index")),
    threads: workflow.cores
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "STAR", "index.time"),
    conda:
        pjoin(ENVS, "star.yaml")
    shell:
        """
        /usr/bin/time -vo {log} STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 149 --genomeSAindexNbases 9
        """


rule STAR_map:
    input:
        index=pjoin(ODIR, "{sample}", "{annov}", "STAR-index"),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    output:
        bam=pjoin(
            ODIR,
            "{sample}",
            "{annov}",
            "STAR",
            "sample_{x}.Aligned.sortedByCoord.out.bam",
        ),
    params:
        oprefix=pjoin(ODIR, "{sample}", "{annov}", "STAR", "sample_{x}."),
    threads: workflow.cores
    log:
        pjoin(ODIR, "{sample}", "bench", "{annov}", "STAR", "map.{x}.time"),
    conda:
        pjoin(ENVS, "star.yaml")
    shell:
        """
        /usr/bin/time -vo {log} STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.oprefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --limitBAMsortRAM 53687091200
        samtools index {output.bam}
        """
