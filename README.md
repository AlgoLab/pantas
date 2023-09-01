# pantas2

Dependencies:
* vg (please compile and use [my branch](https://github.com/ldenti/vg/tree/plain-altid))
* python3-pysam python3-rich python3-biopython
* gffread
* snakemake

``` sh
# Construct, index, and annotate extended splicing graphs
snakemake -s index.smk -c {threads} --config fa=/path/to/reference.fa gtf=/path/to/gene/annotation.gtf vcf=/path/to/phased/variants.vcf.gz odir=/path/to/output/directory
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
