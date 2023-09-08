import sys
from os.path import join as pjoin
from os.path import isfile


FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]
ODIR = config["odir"]

if not isfile(FA + ".fai"):
    print("\n\nInput reference not indexed, please index with samtools faidx\n\n")
    sys.exit(1)

chroms = []
for line in open(FA + ".fai"):
    chroms.append(line.split("\t")[0])

tprefix = ""
for line in open(GTF):
    line = line.strip("\n").split("\t")
    t, info = line[2], line[-1]
    if t in ["mRNA", "transcript"]:
        info = info.split(";")
        for i in info:
            i = i.strip(" ")
            i0, i1 = i.split(" ")
            if i0 == "transcript_id":
                tprefix = i1[1:5]
                break
        break


rule run:
    input:
        pjoin(ODIR, "spliced-pangenome.gcsa"),
        pjoin(ODIR, "spliced-pangenome.gcsa.lcp"),
        pjoin(ODIR, "spliced-pangenome.dist"),
        pjoin(ODIR, "spliced-pangenome.annotated.gfa"),


rule extract_chrom_fa:
    input:
        FA,
    output:
        pjoin(ODIR, "chroms", "{c}.fa"),
    shell:
        """
        samtools faidx {input} {wildcards.c} > {output}
        samtools faidx {output}
        """


rule extract_chrom_gtf:
    input:
        GTF,
    output:
        pjoin(ODIR, "chroms", "{c}.gtf"),
    shell:
        """
        grep -P "^{wildcards.c}\t" {input} > {output}
        """


rule extract_chrom_vcf:
    input:
        VCF,
    output:
        pjoin(ODIR, "chroms", "{c}.vcf.gz"),
    shell:
        """
        bcftools view -Oz {input} {wildcards.c} > {output}
        tabix -p vcf {output}
        """


rule construct:
    input:
        fa=pjoin(ODIR, "chroms", "{c}.fa"),
        vcf=pjoin(ODIR, "chroms", "{c}.vcf.gz"),
    output:
        pg=pjoin(ODIR, "chroms", "{c}", "pangenome.pg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "{c}", "1-construct.txt")
    threads: workflow.cores
    shell:
        """
        vg construct --progress --threads {threads} --reference {input.fa} --vcf {input.vcf} > {output.pg}
        """


rule rna:
    input:
        pg=rules.construct.output.pg,
        gtf=pjoin(ODIR, "chroms", "{c}.gtf"),
    output:
        pg=pjoin(ODIR, "chroms", "{c}", "spliced-pangenome.pg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "{c}", "2-rna.txt")
    threads: workflow.cores
    shell:
        """
        vg rna --progress --threads {threads} --add-ref-paths --transcripts {input.gtf} {input.pg} > {output.pg}
        """


# Prune complex regions and remove alt paths (those starting with _alt_ see https://github.com/vgteam/vg/blob/bcd57125c236782c3f964db10fa581c523ae8e1f/src/path.cpp#L10)
rule prune:
    input:
        pg=rules.rna.output.pg,
    output:
        pg=pjoin(ODIR, "chroms", "{c}", "spliced-pangenome.pruned.pg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "{c}", "3-prune.txt")
    threads: workflow.cores
    shell:
        """
        vg prune --progress --threads {threads} --restore-paths {input.pg} > {output.pg}
        """


# Here the assumption is: thanks to --restore-paths, we have the reference/transcipt paths in the graph, but not the P line
rule reintroduce_paths:
    input:
        gfa=pjoin(ODIR, "chroms", "{c}", "spliced-pangenome.gfa"),
        pgfa=pjoin(ODIR, "chroms", "{c}", "spliced-pangenome.pruned.gfa"),
    output:
        pg=pjoin(ODIR, "chroms", "{c}", "spliced-pangenome.pruned.wpaths.pg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "{c}", "4-reintroducepaths.txt")
    threads: workflow.cores
    shell:
        """
        python3 scripts/reintroduce_paths.py {input.gfa} {input.pgfa} | vg convert --gfa-in --packed-out - > {output.pg}
        """


rule combine:
    input:
        expand(
            pjoin(ODIR, "chroms", "{c}", "spliced-pangenome.pruned.wpaths.pg"),
            c=chroms,
        ),
    output:
        pg=pjoin(ODIR, "spliced-pangenome.wpaths.pg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "5-combine.txt")
    threads: 1
    shell:
        """
        vg combine {input} > {output.pg}
        """


# Keep only reference path as P line (to work properly mpmap needs the reference path **only**)
rule remove_transcripts:
    input:
        pg=pjoin(ODIR, "spliced-pangenome.wpaths.pg"),
    output:
        xg=pjoin(ODIR, "spliced-pangenome.xg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "6-removetranscripts.txt")
    threads: workflow.cores
    shell:
        """
        vg view --threads {threads} {input.pg} | grep -v -P "_R1\t" | vg convert --threads {threads} --gfa-in -x - > {output.xg}
        """


rule index:
    input:
        xg=rules.remove_transcripts.output.xg,
    output:
        gcsa2=pjoin(ODIR, "spliced-pangenome.gcsa"),
        lcp=pjoin(ODIR, "spliced-pangenome.gcsa.lcp"),
        dist=pjoin(ODIR, "spliced-pangenome.dist"),
    params:
        tmpd=ODIR,
    benchmark:
        pjoin(ODIR, "benchmarks", "7-index.txt")
    threads: workflow.cores
    shell:
        """
        vg index --progress --threads {threads} --temp-dir {params.tmpd} --gcsa-out {output.gcsa2} --dist-name {output.dist} {input.xg}
        """


rule annotate:
    input:
        fa=FA,
        gtf=GTF,
        gfa=pjoin(ODIR, "spliced-pangenome.wpaths.gfa"),
    output:
        fa=pjoin(ODIR, "genes.spliced.fa"),
        gfa=pjoin(ODIR, "spliced-pangenome.annotated.gfa"),
    benchmark:
        pjoin(ODIR, "benchmarks", "7-annotate.txt")
    shell:
        """
        gffread -g {input.fa} {input.gtf} -w {output.fa} -W
        python3 scripts/add_junctions.py {input.gfa} {output.fa} > {output.gfa}
        """


rule convert:
    input:
        "{x}.pg",
    output:
        "{x}.gfa",
    shell:
        """
        vg view {input} > {output}
        """
