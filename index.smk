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
        pjoin(ODIR, "index.xg"),
        pjoin(ODIR, "spliced-pangenome.pruned.whaps.annotated.gfa"),


rule construct:
    input:
        fa=FA,
        vcf=VCF,
    output:
        vg=pjoin(ODIR, "pangenome.vg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "1-construct.txt")
    shell:
        """
        vg construct --alt-paths-plain --flat-alts --no-trim-indels --progress --reference {input.fa} --vcf {input.vcf} > {output.vg}
        """


rule rna:
    input:
        vg=rules.construct.output.vg,
        gtf=GTF,
    output:
        vg=pjoin(ODIR, "spliced-pangenome.vg"),
    benchmark:
        pjoin(ODIR, "benchmarks", "2-rna.txt")
    threads: workflow.cores
    shell:
        """
        vg rna --progress --threads {threads} --add-ref-paths --transcripts {input.gtf} {input.vg} > {output.vg}
        """


rule sort_and_convert:
    input:
        vg=rules.rna.output.vg,
    output:
        gfa=pjoin(ODIR, "spliced-pangenome.sorted.gfa"),
    benchmark:
        pjoin(ODIR, "benchmarks", "3-sort.txt")
    shell:
        """
        vg ids --sort {input.vg} | vg view - > {output.gfa}
        """


# FIXME: we need -w 1 or more (depending on multi-allelic sites). Still need to understand why. This can be the reason:
# If first base of a transcript is a bubble, the higher id will be the reference allele.
# All smaller id will be removed from the graph (if w=0)
# so we lose the alternate nodes.
rule prune:
    input:
        gfa=rules.sort_and_convert.output.gfa,
    output:
        gfa=pjoin(ODIR, "spliced-pangenome.pruned.gfa"),
    benchmark:
        pjoin(ODIR, "benchmarks", "4-prune.txt")
    shell:
        """
        python3 scripts/prune_gfa.py -w 3 -t {tprefix} {input.gfa} > {output.gfa}
        """


rule add_haplo:
    input:
        gfa=rules.prune.output.gfa,
        vcf=VCF,
    output:
        gfa=pjoin(ODIR, "spliced-pangenome.pruned.whaps.gfa"),
    benchmark:
        pjoin(ODIR, "benchmarks", "5-addhaplotypes.txt")
    shell:
        """
        python3 scripts/add_haplotypes.py -t {tprefix} {input.gfa} {input.vcf} > {output.gfa}
        """


rule index:
    input:
        gfa=rules.add_haplo.output.gfa,
    output:
        xg=pjoin(ODIR, "index.xg"),
        gcsa2=pjoin(ODIR, "index.gcsa"),
        lcp=pjoin(ODIR, "index.gcsa.lcp"),
    benchmark:
        pjoin(ODIR, "benchmarks", "6-index.txt")
    threads: workflow.cores
    shell:
        """
        vg index --progress --threads {threads} --xg-name {output.xg} --xg-alts --gcsa-out {output.gcsa2} {input.gfa}
        """


rule annotate:
    input:
        fa=FA,
        gtf=GTF,
        gfa=rules.add_haplo.output.gfa,
    output:
        fa=pjoin(ODIR, "genes.spliced.fa"),
        gfa=pjoin(ODIR, "spliced-pangenome.pruned.whaps.annotated.gfa"),
    benchmark:
        pjoin(ODIR, "benchmarks", "7-annotate.txt")
    shell:
        """
        gffread -g {input.fa} {input.gtf} -w {output.fa} -W
        python3 scripts/add_junctions.py {input.gfa} {output.fa} > {output.gfa}
        """
