# pantas

This repository contains the codebase of pantas, a pangenomic approach for performing differential AS events quantification across RNA-Seq conditions. pantas is based on the notion of annotated spliced pangenomes, which are [spliced pangenomes](https://doi.org/10.1038/s41592-022-01731-9) augmented with additional information needed for AS events inference and quantification.

Alongside pantas, we provide a set of utilities to build and index a (annotated) spliced pangenome. We devised [two alternative construction pipelines](https://github.com/AlgoLab/pantas/tree/main#input-preparation), one for full genome analysis and one for the analysis of a panel of genes of interest (reduced).


## Installation
``` sh
git clone https://github.com/AlgoLab/pantas.git

# Install dependencies
mamba create -c bioconda -c conda-forge pantas python=3.10 biopython gffutils intervaltree bcftools samtools gffread vg=1.50.1 snakemake-minimal
mamba activate pantas
```

## pantas pipeline
``` sh
# Augment the annotated spliced pangenome with alignment information (run this for each replicate)
python3 ./scripts/alignments_augmentation_from_gaf.py [condition1-rep1.gaf] [spliced-pangenome.annotated.gfa] > [condition1-rep1.gfa]

# Call events from each graph
python3 ./scripts/call.py [sample.gfa] [annotation.gtf] > [sample-events.csv]

# Quantify events across conditions
python3 ./scripts/quantify3.py -c1 condition1-rep1.csv condition1-rep2.csv condition1-rep3.csv \
                               -c2 condition2-rep1.csv condition2-rep2.csv condition2-rep3.csv > [quantification.csv]
```

### Event calling
The `call.py` script provides several arguments that can be used to tweak the event calling:
```
  --rc RC                        Minimum read count (default: 3)
  --rca RCA                      Minimum read count for annotated events (default: -1)
  --novel                        Call novel events (default: False)
  --no-annotated                 Do not call known annotated events (default: False)
  --events EVENTS [EVENTS ...]   Events to call (default: [ES, SS, IR])
```

To call events from a reduced spliced pangenome, it is necessary to use the `--rp` argument to provide the list of reference paths in the reduced graph (this file is created during the graph construction/indexing):
``` sh
python3 ./scripts/call.py --rp [spliced-pangenome.refpath] [sample.gfa] [annotation.gtf] > [sample-events.csv]
```

### Input preparation
The input of pantas are: an annotated spliced pangenome and the replicates aligned to this graph.

To build and index an annotated spliced pangenome, we provide a snakemake pipeline (`index.smk`): 
``` sh
snakemake -s index.smk -c4 --config fa=/path/to/reference.fa gtf=/path/to/annotation.gtf vcf=/path/to/variants.vcf.gz wd=/path/to/out/dir
```
The annotated spliced pangenome and the index are stored in the `wd` directory:
``` sh
# Annotated spliced pangenome in GFA format:
spliced-pangenome.annotated.gfa
# Compressed graph:
spliced-pangenome.xg
# Index:
spliced-pangenome.dist         
spliced-pangenome.gcsa
spliced-pangenome.gcsa.lcp
```

To build/index a **reduced** annotated spliced pangenomes, i.e., a graph representing a panel of genes of interest:
``` sh
snakemake -s index-reduced.smk -c4 --config fa=/path/to/reference.fa gtf=/path/to/panel.gtf vcf=/path/to/variants.vcf.gz wd=/path/to/out/dir
```
The reduced annotated spliced pangenome and the index are stored in the `wd` directory:
```
# Annotated spliced pangenome in GFA format:
spliced-pangenes.annotated.gfa
# Compressed graph:
spliced-pangenes.xg
# Index:
spliced-pangenes.dist         
spliced-pangenes.gcsa
spliced-pangenes.gcsa.lcp
# Reduced reference paths:
spliced-pangenes.refpath
```

**Note:** using reduced annotated spliced pangenomes (when possible) is recommended since it hugely improves running times and RAM usage.

To map each replicate to the annotated spliced pangenome, we suggest to use `vg mpmap`:
``` sh
vg mpmap -x [spliced-pangenome.xg] -g [spliced-pangenome.gcsa] -d [spliced-pangenome.dist] -f [sample_1.fq] -f [sample_2.fq] -F GAF > [sample.gaf]
```

## Example
The `example` subdirectory contains example data that can be used to test pantas:
``` sh
# Prepare the graph
snakemake -s index.smk -c4 --config fa=example/4.fa gtf=example/4.gtf vcf=example/4.vcf.gz wd=example/pantas-index

# Align the RNA-Seq sample to the graph
vg mpmap -x example/pantas-index/spliced-pangenome.xg \
         -g example/pantas-index/spliced-pangenome.gcsa \
	 -d example/pantas-index/spliced-pangenome.dist \
	 -f example/reads_1.fq -f example/reads_2.fq -F GAF > example/reads.gaf

# Augment the annotated spliced pangenome with alignment information
python3 ./scripts/alignments_augmentation_from_gaf.py example/reads.gaf example/pantas-index/spliced-pangenome.annotated.gfa > example/reads.gfa

# Call all annotated events with minimum support 0 (since example RNA-Seq sample is very small)
# Note that using --rca 0 is equivalent to extract all events from the graph
python3 ./scripts/call.py --rca 0 example/reads.gfa example/4.gtf > example/reads.events.csv

# Quantify the events across the two conditions (an an example here we are using the same file twice)
python3 ./scripts/quantify3.py -c1 example/reads.events.csv -c2 example/reads.events.csv > example/quant.csv
```

## Experiments
Experimental evaluation scripts can be found in the `./exps` subdirectory of this repository. We provide three snakemake pipelines which also contain more information on how to use pantas.
* `./exps/dm-sim/` is the evaluation on simulated data from Drosophila Melanogaster (here we used an annotated spliced pangenome)
* `./exps/dm-sim/` is the evaluation on real data from Drosophila Melanogaster using (here we used an annotated spliced pangenome)
* `./exps/homo-real/` is the evaluation on real data from human (here we used a **reduced** annotated spliced pangenome)


## Authors
pantas is developed by [Simone Ciccolella](https://github.com/sciccolella), [Davide Cozzi](https://github.com/dlcgold), and [Luca Denti](https://github.com/ldenti).

For inquiries on this software please open an [issue](https://github.com/algolab/pantas/issues).
