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
        gfa=pjoin(WD, "pantranscriptome-annotated.gfa"),
        xg=pjoin(WD, "pantranscriptome.xg"),


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
        pg=pjoin(WD, "chroms", "{c}", "pangenome.walts.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "1-construct.txt")
    threads: workflow.cores / 2
    shell:
        """
        vg construct --threads {threads} -r {input.fa} -v {input.vcf} --alt-paths | vg convert --packed-out - > {output.pg}
        """


rule rna_1:
    input:
        pg=rules.construct.output.pg,
        gtf=pjoin(WD, "chroms", "{c}.gtf"),
    output:
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenome.walts.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "2-rna.txt")
    threads: workflow.cores / 4
    shell:
        """
        vg rna --threads {threads} --transcripts {input.gtf} {input.pg} > {output.pg}
        """


rule drop_alts:
    input:
        pg=rules.rna_1.output.pg,
    output:
        pg=pjoin(WD, "chroms", "{c}", "spliced-pangenome.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "3-dropalts.txt")
    threads: 1
    shell:
        """
        vg paths --drop-paths --variant-paths -x {input.pg} | vg convert --packed-out - > {output.pg}
        """


rule gbwt_reference:
    input:
        pg=rules.drop_alts.output.pg,
    output:
        rgbwt=pjoin(WD, "chroms", "{c}", "reference.gbwt"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "4-refgbwt.txt")
    threads: 1
    shell:
        """
        vg gbwt --index-paths --num-jobs {threads} --xg-name {input.pg} --output {output.rgbwt}
        """


rule gbwt:
    input:
        vcf=pjoin(WD, "chroms", "{c}.vcf.gz"),
        pg=rules.rna_1.output.pg,
        rgbwt=rules.gbwt_reference.output.rgbwt,
    output:
        gbwt=pjoin(WD, "chroms", "{c}", "samples.gbwt"),
    params:
        hgbwt=pjoin(WD, "chroms", "{c}", "haplotypes.gbwt"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "5-gbwt.txt")
    threads: workflow.cores
    shell:
        """
        if [[ $(bcftools view -H -G {input.vcf} | head) ]]; then
            vg gbwt --discard-overlaps --num-jobs {threads} --preset 1000gp --vcf-input {input.vcf} --xg-name {input.pg} --output {params.hgbwt}
            vg gbwt --merge -o {output.gbwt} {params.hgbwt} {input.rgbwt}
        else
            cp {input.rgbwt} {output.gbwt}
        fi
        """


rule rna_2:
    input:
        gtf=pjoin(WD, "chroms", "{c}.gtf"),
        pg=rules.drop_alts.output.pg,
        gbwt=rules.gbwt.output.gbwt,
    output:
        txt=pjoin(WD, "chroms", "{c}", "pantranscriptome.info"),
        tgbwt=pjoin(WD, "chroms", "{c}", "pantranscriptome.gbwt"),
        pg=pjoin(WD, "chroms", "{c}", "pantranscriptome.pg"),
    log:
        pjoin(WD, "chroms", "{c}", "pantranscriptome.log"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "6-rna.txt")
    threads: workflow.cores / 2
    shell:
        """
        # --add-ref-paths and --add-hap-paths are needed if we want to prune but keep transcripts
        # but --add-hap-paths makes the indexing too memory intensive (>256GB)
        # --add-ref-paths seems to work on drosophila. If we remove this, we may end up with no reference transcripts and this may break something
        vg rna --progress --add-ref-paths --threads {threads} --haplotypes {input.gbwt} --transcripts {input.gtf} --write-info {output.txt} --write-gbwt {output.tgbwt} {input.pg} > {output.pg} 2> {log}
        """


rule prune:
    input:
        pg=rules.rna_2.output.pg,
    output:
        pg=pjoin(WD, "chroms", "{c}", "pantranscriptome-pruned.pg"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "7-prune.txt")
    log:
        pjoin(WD, "benchmarks", "{c}", "prune.log"),
    threads: workflow.cores / 2
    shell:
        """
        vg prune --progress --threads {threads} --restore-paths {input.pg} > {output.pg} 2> {log}
        """


rule annotate:
    input:
        pg=rules.prune.output.pg,
        txt=rules.rna_2.output.txt,
        gbwt=rules.gbwt.output.gbwt,
        tgbwt=rules.rna_2.output.tgbwt,
    output:
        gfa=pjoin(WD, "chroms", "{c}", "pantranscriptome-annotated.gfa"),
    log:
        pjoin(WD, "chroms", "{c}", "annotation.log"),
    benchmark:
        pjoin(WD, "benchmarks", "{c}", "8-annotate.txt")
    threads: workflow.cores / 2
    shell:
        """
        ./annotate {input.pg} {input.txt} {input.gbwt} {input.tgbwt} > {output.gfa} 2> {log}
        """


rule combine:
    input:
        expand(
            pjoin(WD, "chroms", "{c}", "pantranscriptome-annotated.gfa"),
            c=chroms,
        ),
    output:
        gfa=pjoin(WD, "pantranscriptome-annotated.gfa"),
    benchmark:
        pjoin(WD, "benchmarks", "combine.txt")
    threads: 1
    shell:
        """
        python3 combine.py {input} > {output.gfa}
        """


rule convert:
    input:
        gfa=rules.combine.output.gfa,
    output:
        xg=pjoin(WD, "pantranscriptome.xg"),
    benchmark:
        pjoin(WD, "benchmarks", "convert.txt")
    threads: workflow.cores
    shell:
        """
        vg convert --threads {threads} --gfa-in --xg-out {input.gfa} > {output.xg}
        """
