import sys
from os.path import join as pjoin
from os.path import isfile

FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]
WD = config["wd"]

if not isfile(FA + ".fai"):
    print("\n\nInput reference not indexed, please index with samtools faidx\n\n")
    sys.exit(1)

chroms = []
for line in open(FA + ".fai"):
    chroms.append(line.split("\t")[0])

tprefix = ""
for line in open(GTF):
    if line.startswith("#"):
        continue
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
        pjoin(WD, "spliced-pangenes.gcsa"),
        pjoin(WD, "spliced-pangenes.gcsa.lcp"),
        pjoin(WD, "spliced-pangenes.dist"),
        pjoin(WD, "spliced-pangenes.refpath"),
        pjoin(WD, "spliced-pangenes.annotated.gfa"),


rule extract_chrom_fa:
    input:
        FA,
    output:
        pjoin(WD, "chroms", "{c}.fa"),
    shell:
        """
        samtools faidx {input} {wildcards.c} > {output}
        samtools faidx {output}
        """


rule extract_chrom_gtf:
    input:
        GTF,
    output:
        pjoin(WD, "chroms", "{c}.gtf"),
    shell:
        """
        grep -P "^{wildcards.c}\\t" {input} > {output}
        """


rule extract_chrom_vcf:
    input:
        VCF,
    output:
        pjoin(WD, "chroms", "{c}.vcf.gz"),
    shell:
        """
        bcftools view -Oz {input} {wildcards.c} > {output}
        tabix -p vcf {output}
        """


rule construct:
    input:
        fa=pjoin(WD, "chroms", "{c}.fa"),
        vcf=pjoin(WD, "chroms", "{c}.vcf.gz"),
    output:
        pg=pjoin(WD, "chroms", "{c}", "pangenome.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "1-construct.txt")
    threads: workflow.cores
    shell:
        """
        vg construct --progress --threads {threads} --reference {input.fa} --vcf {input.vcf} > {output.pg}
        """


rule rna:
    input:
        pg=rules.construct.output.pg,
        gtf=pjoin(WD, "chroms", "{c}.gtf"),
    output:
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenome.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "2-rna.txt")
    threads: workflow.cores
    shell:
        """
        vg rna --progress --threads {threads} --add-ref-paths --transcripts {input.gtf} {input.pg} > {output.pg}
        """


rule tsort:
    input:
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenome.pg"),
    output:
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenome.tsorted.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "2a-tsort.txt")
    threads: 1
    shell:
        """
        vg ids -s {input.pg} > {output.pg}
        """


# Prune complex regions and remove alt paths (those starting with _alt_ see https://github.com/vgteam/vg/blob/bcd57125c236782c3f964db10fa581c523ae8e1f/src/path.cpp#L10)
rule prune:
    input:
        pg=rules.tsort.output.pg,
    output:
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenome.pruned.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "3-prune.txt")
    threads: workflow.cores
    shell:
        """
        vg prune --progress --threads {threads} --restore-paths {input.pg} > {output.pg}
        """


# Here the assumption is: thanks to --restore-paths, we have the reference/transcript paths in the graph, but not the P line
rule reintroduce_paths:
    input:
        gfa=pjoin(WD, "chroms", "{c}", "spliced-pangenome.tsorted.gfa"),
        pgfa=pjoin(WD, "chroms", "{c}", "spliced-pangenome.pruned.gfa"),
    output:
        gfa=pjoin(WD, "chroms", "{c}", "spliced-pangenome.pruned.wpaths.gfa"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "4-reintroducepaths.txt")
    threads: workflow.cores
    shell:
        """
        python3 scripts/reintroduce_paths.py {input.gfa} {input.pgfa} > {output.gfa}
        """

rule reduce_graph:
    input:
        gfa=rules.reintroduce_paths.output.gfa,
    output:
        gfa=pjoin(WD, "chroms", "{c}", "spliced-pangenes.gfa"),
        gfa_rf=pjoin(WD, "chroms", "{c}", "spliced-pangenes.gfa.refpath"),
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenes.pg"),
    shell:
        """
        python3 scripts/reduce.py {input.gfa} {output.gfa}
        vg convert --gfa-in {output.gfa} --packed-out > {output.pg}
        """


rule combine:
    input:
        expand(
            pjoin(WD, "chroms", "{c}", "spliced-pangenes.pg"),
            c=chroms,
        ),
    output:
        pg=pjoin(WD, "spliced-pangenes.wpaths.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "5-combine.txt")
    threads: 1
    shell:
        """
        vg combine {input} > {output.pg}
        """

rule combine_refpaths:
    input:
        expand(
            pjoin(WD, "chroms", "{c}", "spliced-pangenes.gfa.refpath"),
            c=chroms,
        ),
    output:
        gfa_rp=pjoin(WD, "spliced-pangenes.refpath"),
    shell:
        """
        cat {input} > {output.gfa_rp}
        """


# Keep only reference path as P line (to work properly mpmap needs the reference paths **only**)
rule remove_transcripts:
    input:
        pg=rules.combine.output.pg,
    output:
        xg=pjoin(WD, "spliced-pangenes.xg"),
    benchmark:
        pjoin(WD, "benchmarks", "6-removetranscripts.txt")
    threads: workflow.cores
    shell:
        """
        vg view --threads {threads} {input.pg} | grep -v -P "_R1\\t" | vg convert --threads {threads} --gfa-in -x - > {output.xg}
        """


rule index:
    input:
        xg=rules.remove_transcripts.output.xg,
    output:
        gcsa2=pjoin(WD, "spliced-pangenes.gcsa"),
        lcp=pjoin(WD, "spliced-pangenes.gcsa.lcp"),
        dist=pjoin(WD, "spliced-pangenes.dist"),
    params:
        tmpd=WD,
    benchmark:
        pjoin(WD, "benchmarks", "7-index.txt")
    threads: workflow.cores
    shell:
        """
        vg index --progress --threads {threads} --temp-dir {params.tmpd} --gcsa-out {output.gcsa2} --dist-name {output.dist} {input.xg}
        """


rule annotate:
    input:
        fa=FA,
        gtf=GTF,
        gfa=pjoin(WD, "spliced-pangenes.wpaths.gfa"),
    output:
        fa=pjoin(WD, "genes.spliced.fa"),
        gfa=pjoin(WD, "spliced-pangenes.annotated.gfa"),
    benchmark:
        pjoin(WD, "benchmarks", "7-annotate.txt")
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
