rule STAR_index_anno:
    input:
        fa=FA,
        gtf=GTF,
    output:
        index=directory(pjoin(ODIR, "STAR-index")),
    threads: workflow.cores
    conda:
        "../envs/star.yaml"
    log:
        pjoin(ODIR, "bench", "STAR", "index.time"),
    shell:
        """
        /usr/bin/time -vo {log} STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 149 --genomeSAindexNbases 9
        """


rule STAR_map:
    input:
        index=pjoin(ODIR, "STAR-index"),
        fq1=lambda wildcards: FQs[wildcards.sample][0],
        fq2=lambda wildcards: FQs[wildcards.sample][1],
    output:
        bam=pjoin(ODIR, "{sample}.STAR/Aligned.sortedByCoord.out.bam"),
    params:
        oprefix=pjoin(ODIR, "{sample}.STAR/"),
    threads: workflow.cores
    conda:
        "../envs/star.yaml"
    log:
        pjoin(ODIR, "bench", "STAR", "map.{sample}.time"),
    shell:
        """
        /usr/bin/time -vo {log} STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.oprefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --limitBAMsortRAM 53687091200
        samtools index {output.bam}
        """


rule rmats_c1:
    input:
        expand(
            pjoin(ODIR, "{sample}.STAR/Aligned.sortedByCoord.out.bam"),
            sample=C1.keys(),
        ),
    output:
        pjoin(ODIR, "rMATS", "c1-bams.txt"),
    shell:
        """
        echo {input} | tr " " "," > {output}
        """


rule rmats_c2:
    input:
        expand(
            pjoin(ODIR, "{sample}.STAR/Aligned.sortedByCoord.out.bam"),
            sample=C2.keys(),
        ),
    output:
        pjoin(ODIR, "rMATS", "c2-bams.txt"),
    shell:
        """
        echo {input} | tr " " "," > {output}
        """


rule rmats:
    input:
        gtf=GTF,
        c1txt=pjoin(ODIR, "rMATS", "c1-bams.txt"),
        c2txt=pjoin(ODIR, "rMATS", "c2-bams.txt"),
    output:
        summary=pjoin(ODIR, "rMATS", "summary.txt"),
    params:
        tmpd=pjoin(ODIR, "rMATS-tmp"),
        outd=pjoin(ODIR, "rMATS"),
    threads: workflow.cores
    conda:
        "../envs/rmats.yaml"
    log:
        pjoin(ODIR, "bench", "rmats.time"),
    shell:
        """
        /usr/bin/time -vo {log} rmats.py --gtf {input.gtf} --b1 {input.c1txt} --b2 {input.c2txt} --od {params.outd} --tmp {params.tmpd} --readLength {L} --nthread {threads} -t paired
        """
