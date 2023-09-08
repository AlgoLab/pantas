# pantas2

Dependencies:
* vg
* python3-pysam python3-rich python3-biopython
* gffread
* snakemake

``` sh
snakemake -s index.smk -c4 --config fa=example/4.fa gtf=example/4.gtf vcf=example/4.vcf.gz odir=example/OUT

# TODO: add haplotyoes somehow?

# Map the reads
vg mpmap -x example/OUT/spliced-pangenome.xg -g example/OUT/spliced-pangenome.gcsa -d example/OUT/spliced-pangenome.dist -f example/reads_1.fq -f example/reads_2.fq -F GAF > example/reads.gaf

# Weight the graph
python3 ./scripts/alignments_augmentation.py example/reads.gaf x example/OUT/spliced-pangenome.annotated.gfa > example/spliced-pangenome.weighted.gfa
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
