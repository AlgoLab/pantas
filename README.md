# pantas2

Dependencies:
* vg (please compile and use [my branch](https://github.com/ldenti/vg/tree/plain-altid))
* python3-pysam python3-rich python3-biopython
* gffread
* snakemake

``` sh
# Construct and index
vg autoindex --workflow mpmap --prefix example/4-index --ref-fasta example/4.fa --vcf example/4.vcf.gz --tx-gff example/4.gtf --tmp-dir . --threads 4 --verbosity 2
# Extract graph with transcripts as paths
vg convert --packed-out example/4-index.spliced.xg | vg rna --progress --threads 16 --add-ref-paths --transcripts example/4.gtf - | vg view - > example/4-graph.gfa
# Annotate graph
gffread -g example/4.fa example/4.gtf -w example/4.genes-spliced.fa -W
python3 ./scripts/add_junctions.py example/4-graph.gfa example/4.genes-spliced.fa > example/4-graph.anno.gfa

# Align sample
vg mpmap -x example/4-index.spliced.xg -g example/4-index.spliced.gcsa -d example/4-index.spliced.dist -f example/reads_1.fq -f example/reads_2.fq -F GAMP > example/reads.gamp
vg view -K -j example/reads.gamp > example/reads.json
python3 ./scripts/alignments_augmentation.py example/reads.json x example/4-graph.anno.gfa > example/4-graph.anno.weighted.gfa
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
