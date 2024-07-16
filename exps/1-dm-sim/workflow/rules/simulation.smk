rule apply_variants:
    input:
        fa=FA,
        vcf=VCF,
    output:
        fa=pjoin(ODIR, "{sample}", "ref-wvars.fa"),
    shell:
        """
        bcftools consensus -f {input.fa} -H 1 -s {wildcards.sample} {input.vcf} > {output.fa}
        samtools faidx {output.fa}
        """


rule extract_chrom:
    input:
        fa=rules.apply_variants.output.fa,
    output:
        fa=pjoin(ODIR, "{sample}", "asim-input", "{c}.fa"),
    shell:
        """
        samtools faidx {input.fa} {wildcards.c} > {output.fa}
        """


rule link_annotation:
    input:
        GTF,
    output:
        pjoin(ODIR, "{sample}", "asim-input", "annotation.gtf"),
    shell:
        """
        ln {input} {output}
        """


rule simulate:
    input:
        expand(pjoin(ODIR, "{{sample}}", "asim-input", "{c}.fa"), c=chroms),
        pjoin(ODIR, "{sample}", "asim-input", "annotation.gtf"),
    output:
        pjoin(ODIR, "{sample}", "asim-output", "sample_01_1.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "sample_01_2.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "sample_02_1.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "sample_02_2.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "exon_junction_coverage.tsv"),
        pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf"),
        pjoin(ODIR, "{sample}", "asim-output", "event_annotation.tsv"),
    params:
        input_dir=pjoin(ODIR, "{sample}", "asim-input"),
        output_dir=pjoin(ODIR, "{sample}", "asim-output"),
    threads: workflow.cores
    shell:
        """
        Rscript scripts/asimulator.R {params.input_dir} {params.output_dir} {L} {N} {seed} {threads}
        """


rule get_novel:
    input:
        gtf=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf"),
    output:
        gtf=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf"),
    shell:
        """
        grep -v 'template "FALSE"' {input.gtf} > {output.gtf}
        """


# --------------
# --- TRUTH  ---
# --------------
rule clean_samples:
    input:
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.fastq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.fastq"),
    output:
        pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_2.clean.fq"),
    shell:
        """
        python3 scripts/filter_reads.py {input.fq1} {input.fq2}
        """


rule get_rc:
    input:
        tsv1=pjoin(ODIR, "{sample}", "asim-output", "exon_junction_coverage.tsv"),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "sample_0{x}_1.clean.fq"),
        tsv2=pjoin(ODIR, "{sample}", "asim-output", "event_annotation.tsv"),
    output:
        csv=pjoin(ODIR, "{sample}", "asim-output", "read-counts.{x}.csv"),
    shell:
        """
        python3 scripts/simrc.py {input.fq1} {input.tsv1} {input.tsv2} > {output.csv}
        """


rule merge_rc:
    input:
        csv1=pjoin(ODIR, "{sample}", "asim-output", "read-counts.1.csv"),
        csv2=pjoin(ODIR, "{sample}", "asim-output", "read-counts.2.csv"),
    output:
        csv=pjoin(ODIR, "{sample}", "asim-output", "read-counts.csv"),
    shell:
        """
        # CHECKME: assuming ordered input
        paste {input.csv1}  <(cut -f11 -d',' {input.csv2} ) -d',' > {output.csv}
        """


rule get_truth:
    input:
        tsv=pjoin(ODIR, "{sample}", "asim-output", "event_annotation.tsv"),
        csv=pjoin(ODIR, "{sample}", "asim-output", "read-counts.csv"),
    output:
        csv=pjoin(ODIR, "{sample}", "truth.csv"),
    shell:
        """
        python3 scripts/build_truth.py {input.tsv} {input.csv} > {output.csv}
        """
