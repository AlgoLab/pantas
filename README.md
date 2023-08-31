# pantas2

Dependencies:
* vg (please compile and use [my brach](https://github.com/ldenti/vg/tree/plain-altid))
* python3-pysam python3-rich

``` sh
# Filter VCF by genes
grep -P "\tgene\t" chr21.gtf | cut -f1,4,5 > chr21.genes.bed
bcftools view -R chr21.pcg.bed -c 1 -q 0.001 chr21.genes.vcf.gz  | bcftools norm -Oz -f chr21.fa > chr21.pcg.vcf.gz

# Construct and index graph
vg construct --alt-paths-plain --reference chr21.fa --vcf chr21.pcg.vcf.gz --flat-alts --no-trim-indels --progress > pcg.vg
vg rna --progress --threads 4 --add-ref-paths --transcripts chr21.pcg.gtf pcg.vg > pcg.spliced.vg
vg ids --sort pcg.spliced.vg  > pcg.spliced.sorted.vg
vg view pcg.spliced.sorted.vg  > pcg.spliced.sorted.gfa
python3 prune.py -w 0 -t ENST pcg.spliced.sorted.gfa > pcg.spliced.pruned.gfa
python3 add_haplotypes.py -t ENST pcg.spliced.pruned.gfa chr21.pcg.vcf.gz > pcg.spliced.pruned.whaps.gfa
vg index --progress --threads 4 --xg-name pcg.index.xg --xg-alts --gcsa-out pcg.index.gcsa pcg.spliced.pruned.whaps.gfa
```
