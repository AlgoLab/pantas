rule apply_variants:
    input:
        fa=FA,
        vcf=VCF,
    output:
        fa=pjoin(ODIR, "{sample}", "ref-wvars-{h}.fa"),
    shell:
        """
        bcftools consensus -f {input.fa} -H {wildcards.h} -s {wildcards.sample} {input.vcf} > {output.fa}
        samtools faidx {output.fa}
        """


rule extract_chrom:
    input:
        fa=rules.apply_variants.output.fa,
    output:
        fa=pjoin(ODIR, "{sample}", "asim-input", "{h}", "{c}.fa"),
    shell:
        """
        samtools faidx {input.fa} {wildcards.c} > {output.fa}
        """


rule link_annotation:
    input:
        GTF,
    output:
        pjoin(ODIR, "{sample}", "asim-input", "{h}", "annotation.gtf"),
    shell:
        """
        ln {input} {output}
        """


rule simulate:
    input:
        expand(pjoin(ODIR, "{{sample}}", "asim-input", "{{h}}", "{c}.fa"), c=chroms),
        pjoin(ODIR, "{sample}", "asim-input", "{h}", "annotation.gtf"),
    output:
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_01_1.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_01_2.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_02_1.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_02_2.fastq"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "exon_junction_coverage.tsv"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "splicing_variants.gtf"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "event_annotation.tsv"),
    params:
        input_dir=pjoin(ODIR, "{sample}", "asim-input", "{h}"),
        output_dir=pjoin(ODIR, "{sample}", "asim-output", "{h}"),
    threads: workflow.cores
    shell:
        """
        Rscript scripts/asimulator.R {params.input_dir} {params.output_dir} {L} {N} {seed} {threads}
        """


rule get_novel:
    input:
        gtf=pjoin(ODIR, "{sample}", "asim-output", "{h}", "splicing_variants.gtf"),
    output:
        gtf=pjoin(ODIR, "{sample}", "asim-output", "{h}", "splicing_variants_novel.gtf"),
    shell:
        """
        grep -v 'template "FALSE"' {input.gtf} > {output.gtf}
        """


# --------------
# --- TRUTH  ---
# --------------
rule clean_samples:
    input:
        fq1=pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_0{x}_1.fastq"),
        fq2=pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_0{x}_2.fastq"),
    output:
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_0{x}_1.clean.fq"),
        pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_0{x}_2.clean.fq"),
    shell:
        """
        python3 scripts/filter_reads.py {input.fq1} {input.fq2}
        """


rule get_rc:
    input:
        tsv1=pjoin(ODIR, "{sample}", "asim-output", "{h}", "exon_junction_coverage.tsv"),
        fq1=pjoin(ODIR, "{sample}", "asim-output", "{h}", "sample_0{x}_1.clean.fq"),
        tsv2=pjoin(ODIR, "{sample}", "asim-output", "{h}", "event_annotation.tsv"),
    output:
        csv=pjoin(ODIR, "{sample}", "asim-output", "{h}", "read-counts.{x}.csv"),
    shell:
        """
        python3 scripts/simrc.py {input.fq1} {input.tsv1} {input.tsv2} > {output.csv}
        """


rule merge_rc:
    input:
        csv11=pjoin(ODIR, "{sample}", "asim-output", "1", "read-counts.1.csv"),
        csv12=pjoin(ODIR, "{sample}", "asim-output", "1", "read-counts.2.csv"),
        csv21=pjoin(ODIR, "{sample}", "asim-output", "2", "read-counts.1.csv"),
        csv22=pjoin(ODIR, "{sample}", "asim-output", "2", "read-counts.2.csv"),
    output:
        csv=pjoin(ODIR, "{sample}", "asim-output", "read-counts.csv"),
    shell:
        """
        python3 scripts/merge_rc.py {input.csv11} {input.csv21} {input.csv12} {input.csv22} > {output.csv}
        """


rule cat_fq:
    input:
        fq11_1=pjoin(ODIR, "{sample}", "asim-output", "1", "sample_01_1.clean.fq"),
        fq11_2=pjoin(ODIR, "{sample}", "asim-output", "1", "sample_01_2.clean.fq"),
        fq12_1=pjoin(ODIR, "{sample}", "asim-output", "1", "sample_02_1.clean.fq"),
        fq12_2=pjoin(ODIR, "{sample}", "asim-output", "1", "sample_02_2.clean.fq"),
        fq21_1=pjoin(ODIR, "{sample}", "asim-output", "2", "sample_01_1.clean.fq"),
        fq21_2=pjoin(ODIR, "{sample}", "asim-output", "2", "sample_01_2.clean.fq"),
        fq22_1=pjoin(ODIR, "{sample}", "asim-output", "2", "sample_02_1.clean.fq"),
        fq22_2=pjoin(ODIR, "{sample}", "asim-output", "2", "sample_02_2.clean.fq"),
    output:
        fq1_1=pjoin(ODIR, "{sample}", "asim-output", "sample_01_1.clean.fq"),
        fq1_2=pjoin(ODIR, "{sample}", "asim-output", "sample_01_2.clean.fq"),
        fq2_1=pjoin(ODIR, "{sample}", "asim-output", "sample_02_1.clean.fq"),
        fq2_2=pjoin(ODIR, "{sample}", "asim-output", "sample_02_2.clean.fq"),
    shell:
        """
        cat {input.fq11_1} {input.fq21_1} > {output.fq1_1}
        cat {input.fq11_2} {input.fq21_2} > {output.fq1_2}
        cat {input.fq12_1} {input.fq22_1} > {output.fq2_1}
        cat {input.fq12_2} {input.fq22_2} > {output.fq2_2}
        """


rule cp_files_from_1:
    input:
        tsv=pjoin(ODIR, "{sample}", "asim-output", "1", "event_annotation.tsv"),
        gtf1=pjoin(ODIR, "{sample}", "asim-output", "1", "splicing_variants.gtf"),
        gtf2=pjoin(ODIR, "{sample}", "asim-output", "1", "splicing_variants_novel.gtf"),
    output:
        tsv=pjoin(ODIR, "{sample}", "asim-output", "event_annotation.tsv"),
        gtf1=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants.gtf"),
        gtf2=pjoin(ODIR, "{sample}", "asim-output", "splicing_variants_novel.gtf"),
    shell:
        """
        cp {input.tsv} {output.tsv}
        cp {input.gtf1} {output.gtf1}
        cp {input.gtf2} {output.gtf2}
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
