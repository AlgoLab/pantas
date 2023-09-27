import sys
import random
from os.path import join as pjoin
from os.path import isfile
import gzip

# FIXME: everything is hardcoded to 1 replicate and 2 conditions

random.seed(23)

FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]
ODIR = config["odir"]
L = 150  # config["l"]
N = config["n"]

if not isfile(FA + ".fai"):
    print("\n\nInput reference not indexed, please index with samtools faidx\n\n")
    sys.exit(1)

chroms = []
for line in open(FA + ".fai"):
    chroms.append(line.split("\t")[0])

sample = None
with gzip.open(VCF, "rt") as fin:
    for line in fin:
        if line.startswith("#C"):
            samples = line.split("\t")[9:]
            sample = random.choice(samples)
            break
assert sample != None


rule run:
    input:
        pjoin(ODIR, sample, "truth.csv"),
        expand(pjoin(ODIR, sample, "pantas2-anno", "sample_{x}.bam"), x=[1, 2]),
        pjoin(ODIR, sample, "pantas2-anno", "quant.csv"),
        pjoin(ODIR, sample, "rMATS", "summary.txt"),


# ------------------
# --- SIMULATION ---
# ------------------
rule apply_variants:
    input:
        fa=FA,
        vcf=VCF,
    output:
        fa=pjoin(ODIR, f"{sample}.fa"),
    shell:
        """
        bcftools consensus -f {input.fa} -H 1 -s {sample} {input.vcf} > {output.fa}
        samtools faidx {output.fa}
        """


rule extract_chrom:
    input:
        fa=rules.apply_variants.output.fa,
    output:
        fa=pjoin(ODIR, sample, "asim-input", "{c}.fa"),
    shell:
        """
        samtools faidx {input.fa} {wildcards.c} > {output.fa}
        """


rule link_annotation:
    input:
        GTF,
    output:
        pjoin(ODIR, sample, "asim-input", "annotation.gtf"),
    shell:
        """
        ln {input} {output}
        """


rule simulate:
    input:
        expand(pjoin(ODIR, sample, "asim-input", "{c}.fa"), c=chroms),
        pjoin(ODIR, sample, "asim-input", "annotation.gtf"),
    output:
        pjoin(ODIR, sample, "asim-output", "sample_01_1.fastq"),
        pjoin(ODIR, sample, "asim-output", "sample_01_2.fastq"),
        pjoin(ODIR, sample, "asim-output", "sample_02_1.fastq"),
        pjoin(ODIR, sample, "asim-output", "sample_02_2.fastq"),
        pjoin(ODIR, sample, "asim-output", "exon_junction_coverage.tsv"),
        pjoin(ODIR, sample, "asim-output", "splicing_variants.gtf"),
        pjoin(ODIR, sample, "asim-output", "splicing_variants_novel.gtf"),
        pjoin(ODIR, sample, "asim-output", "event_annotation.tsv"),
    params:
        input_dir=pjoin(ODIR, sample, "asim-input"),
        output_dir=pjoin(ODIR, sample, "asim-output"),
    shell:
        """
        Rscript asimulator.R {params.input_dir} {params.output_dir} {L} {N}
        """


# --------------
# --- TRUTH  ---
# --------------
rule clean_samples:
    input:
        fq1=pjoin(ODIR, sample, "asim-output", "sample_0{x}_1.fastq"),
        fq2=pjoin(ODIR, sample, "asim-output", "sample_0{x}_2.fastq"),
    output:
        pjoin(ODIR, sample, "asim-output", "sample_0{x}_1.clean.fq"),
        pjoin(ODIR, sample, "asim-output", "sample_0{x}_2.clean.fq"),
    shell:
        """
        python3 filter_reads.py {input.fq1} {input.fq2}
        """


rule get_rc:
    input:
        tsv1=pjoin(ODIR, sample, "asim-output", "exon_junction_coverage.tsv"),
        fq1=pjoin(ODIR, sample, "asim-output", "sample_0{x}_1.clean.fq"),
        tsv2=pjoin(ODIR, sample, "asim-output", "event_annotation.tsv"),
    output:
        csv=pjoin(ODIR, sample, "asim-output", "read-counts.{x}.csv"),
    shell:
        """
        python3 simrc.py {input.fq1} {input.tsv1} {input.tsv2} > {output.csv}
        """


rule merge_rc:
    input:
        csv1=pjoin(ODIR, sample, "asim-output", "read-counts.1.csv"),
        csv2=pjoin(ODIR, sample, "asim-output", "read-counts.2.csv"),
    output:
        csv=pjoin(ODIR, sample, "asim-output", "read-counts.csv"),
    shell:
        """
        # CHECKME: assuming ordered input
        paste {input.csv1}  <(cut -f11 -d',' {input.csv2} ) -d',' > {output.csv}
        """


rule get_truth:
    input:
        tsv=pjoin(ODIR, sample, "asim-output", "event_annotation.tsv"),
        csv=pjoin(ODIR, sample, "asim-output", "read-counts.csv"),
    output:
        csv=pjoin(ODIR, sample, "truth.csv"),
    shell:
        """
        python3 build_truth.py {input.tsv} {input.csv} > {output.csv}
        """


# --------------
# --- PANTAS ---
# --------------
rule pantas2_index_anno:
    input:
        fa=FA,
        gtf=pjoin(ODIR, sample, "asim-output", "splicing_variants.gtf"),
        vcf=VCF,
    output:
        xg=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.dist"),
        gfa=pjoin(
            ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.annotated.gfa"
        ),
        splicedfa=pjoin(ODIR, sample, "pantas2-anno", "index", "genes.spliced.fa"),
    params:
        wd=pjoin(ODIR, sample, "pantas2-anno", "index"),
    log:
        time=pjoin(ODIR, sample, "pantas2-anno", "index.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} bash ../pantas2.sh index {input.fa} {input.gtf} {input.vcf} {params.wd} {threads}
        """


# FIXME: use wildcard for anno/novel
rule pantas2_mpmap_anno:
    input:
        xg=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.dist"),
        fq1=pjoin(ODIR, sample, "asim-output", "sample_0{x}_1.fastq"),
        fq2=pjoin(ODIR, sample, "asim-output", "sample_0{x}_2.fastq"),
    output:
        gaf=pjoin(ODIR, sample, "pantas2-anno", "sample_{x}.gaf"),
    threads: workflow.cores
    log:
        pjoin(ODIR, sample, "pantas2-anno", "mpmap{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F GAF --threads {threads} > {output.gaf}
        """


rule pantas2_mpmap_bam_anno:
    input:
        xg=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.xg"),
        gcsa=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.gcsa"),
        dist=pjoin(ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.dist"),
        fq1=pjoin(ODIR, sample, "asim-output", "sample_0{x}_1.fastq"),
        fq2=pjoin(ODIR, sample, "asim-output", "sample_0{x}_2.fastq"),
    output:
        bam=pjoin(ODIR, sample, "pantas2-anno", "sample_{x}.bam"),
    threads: workflow.cores
    shell:
        """
        vg mpmap -x {input.xg} -g {input.gcsa} -d {input.dist} -f {input.fq1} -f {input.fq2} -F BAM --threads {threads} | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule pantas_weight_anno:
    input:
        gfa=pjoin(
            ODIR, sample, "pantas2-anno", "index", "spliced-pangenome.annotated.gfa"
        ),
        gaf=pjoin(ODIR, sample, "pantas2-anno", "sample_{x}.gaf"),
    output:
        gfa=pjoin(ODIR, sample, "pantas2-anno", "graph_{x}.gfa"),
    log:
        pjoin(ODIR, sample, "pantas2-anno", "weight{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 ../scripts/alignments_augmentation_from_gaf.py {input.gaf} {input.gfa} > {output.gfa}
        """


rule pantas_call_anno:
    input:
        gfa=pjoin(ODIR, sample, "pantas2-anno", "graph_{x}.gfa"),
        gtf=pjoin(ODIR, sample, "asim-output", "splicing_variants.gtf"),
    output:
        csv=pjoin(ODIR, sample, "pantas2-anno", "events_{x}.csv"),
    log:
        pjoin(ODIR, sample, "pantas2-anno", "call{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 ../scripts/call.py --rc -1 {input.gfa} {input.gtf} > {output.csv}
        """


rule pantas_quant_anno:
    input:
        splicedfa=pjoin(ODIR, sample, "pantas2-anno", "index", "genes.spliced.fa"),
        csv1=pjoin(ODIR, sample, "pantas2-anno", "events_1.csv"),
        csv2=pjoin(ODIR, sample, "pantas2-anno", "events_2.csv"),
    output:
        csv=pjoin(ODIR, sample, "pantas2-anno", "quant.csv"),
    log:
        pjoin(ODIR, sample, "pantas2-anno", "quant.time"),
    shell:
        """
        /usr/bin/time -vo {log} python3 ../scripts/quantify.py {input.splicedfa} {input.csv1} {input.csv2} > {output.csv}
        """


# --------------------
# --- STAR + RMATS ---
# --------------------
rule STAR_index:
    input:
        fa=FA,
        gtf=pjoin(ODIR, sample, "asim-output", "splicing_variants.gtf"),
    output:
        index=directory(pjoin(ODIR, sample, "STAR-index")),
    threads: workflow.cores
    conda:
        "envs/star.yaml"
    log:
        pjoin(ODIR, sample, "bench", "star-index.time"),
    shell:
        """
        /usr/bin/time -vo {log} STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 149 --genomeSAindexNbases 9
        """


rule STAR_map:
    input:
        index=pjoin(ODIR, sample, "STAR-index"),
        fq1=pjoin(ODIR, sample, "asim-output", "sample_0{x}_1.fastq"),
        fq2=pjoin(ODIR, sample, "asim-output", "sample_0{x}_2.fastq"),
    output:
        bam=pjoin(ODIR, sample, "STAR", "sample_{x}.Aligned.sortedByCoord.out.bam"),
    params:
        oprefix=pjoin(ODIR, sample, "STAR", "sample_{x}."),
    threads: workflow.cores
    conda:
        "envs/star.yaml"
    log:
        pjoin(ODIR, sample, "bench", "star-map.{x}.time"),
    shell:
        """
        /usr/bin/time -vo {log} STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix {params.oprefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --limitBAMsortRAM 53687091200
        samtools index {output.bam}
        """


rule rmats_c1:
    input:
        pjoin(ODIR, sample, "STAR", "sample_1.Aligned.sortedByCoord.out.bam"),
    output:
        pjoin(ODIR, sample, "STAR.1.txt"),
    shell:
        """
        echo {input} > {output}
        """


rule rmats_c2:
    input:
        pjoin(ODIR, sample, "STAR", "sample_2.Aligned.sortedByCoord.out.bam"),
    output:
        pjoin(ODIR, sample, "STAR.2.txt"),
    shell:
        """
        echo {input} > {output}
        """


rule rmats:
    input:
        gtf=pjoin(ODIR, sample, "asim-output", "splicing_variants.gtf"),
        c1txt=pjoin(ODIR, sample, "STAR.1.txt"),
        c2txt=pjoin(ODIR, sample, "STAR.2.txt"),
    output:
        outd=directory(pjoin(ODIR, sample, "rMATS")),
        summary=pjoin(ODIR, sample, "rMATS", "summary.txt"),
    params:
        tmpd=pjoin(ODIR, sample, "rMATS-tmp"),
    threads: workflow.cores
    conda:
        "envs/rmats.yaml"
    log:
        pjoin(ODIR, sample, "bench", "rmats.time"),
    shell:
        """
        /usr/bin/time -vo {log} rmats.py --gtf {input.gtf} --b1 {input.c1txt} --b2 {input.c2txt} --od {output.outd} --tmp {params.tmpd} --readLength 150 --nthread {threads} -t paired
        """
