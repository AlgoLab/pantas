# pantas2

Dependencies:
* vg (please compile and use [my branch](https://github.com/ldenti/vg/tree/plain-altid))
* python3-pysam python3-rich python3-biopython
* gffread
* snakemake

``` sh
# Construct and index graph
vg construct --alt-paths-plain --reference chr21.fa --vcf chr21.pcg.vcf.gz --flat-alts --no-trim-indels --progress > pcg.vg
vg rna --progress --threads 4 --add-ref-paths --transcripts chr21.pcg.gtf pcg.vg > pcg.spliced.vg
vg ids --sort pcg.spliced.vg  > pcg.spliced.sorted.vg
vg view pcg.spliced.sorted.vg  > pcg.spliced.sorted.gfa
python3 prune.py -w 0 -t ENST pcg.spliced.sorted.gfa > pcg.spliced.pruned.gfa
python3 add_haplotypes.py -t ENST pcg.spliced.pruned.gfa chr21.pcg.vcf.gz > pcg.spliced.pruned.whaps.gfa
vg index --progress --threads 4 --xg-name pcg.index.xg --xg-alts --gcsa-out pcg.index.gcsa pcg.spliced.pruned.whaps.gfa

# Annotate GFA with Exons, Transcripts and Junctions 
# (these commands apply to the provided test files)
gffread -g test/reference.fa test/gene.gtf -w test/genes.spliced.fa -W
python scripts/add_junctions.py test/graph.spliced.pruned.whaps.gfa test/genes.spliced.fa > test/graph.spliced.pruned.whaps.annotated.gfa
```

## New GFA fields

### `S` Segment line

| Tag 	| Type 	| Description 	|
|---	|---	|---	|
| EX 	| Z 	| Comma (`,`) separated values for each exons the node is part of. Annotated as `[Transcript].[Exon_Number]` 	|

##### Example:
```
S	5	AAA	LN:i:3	EX:Z:Ttest.1
S       577768  ATTTAT  LN:i:6
S       549841  TTCATCTGGTAGTTCTTG      LN:i:18 EX:Z:ENST00000284878.4,ENST00000400166.4,ENST00000400165.4,ENST00000400169.4
```

### `L` Link line

| Tag 	| Type 	| Description 	|
|---	|---	|---	|
| JN 	| Z 	| Comma (`,`) separated values for each junction between exons. Annotated as `[Transcript].[Exon_From].[Exon_To]` 	|

##### Example:
```
L	15	+	16	+	0M	JN:Z:Ttest.2.3
L       548094  +       548095  +       0M
L       580290  +       580291  +       0M      JN:Z:ENST00000284885.6.7,ENST00000422787.7.8,ENST00000474775.3.4
```
