# pantas

This repository contains the codebase of pantas, a pangenomic approach for performing differential AS events quantification across RNA-Seq conditions. pantas is based on the notion of annotated spliced pangenomes, which are [spliced pangenomes](https://doi.org/10.1038/s41592-022-01731-9) augmented with additional information needed for AS events inference and quantification.

Alongside pantas, we provide a set of utilities to build and index a (annotated) spliced pangenome. We devised [two alternative construction pipelines](https://github.com/AlgoLab/pantas/tree/main#input-preparation), one for full genome analysis and one for the analysis of a panel of genes of interest (reduced indexing).


## Installation
``` sh
git clone https://github.com/AlgoLab/pantas.git

# Install dependencies (all available from bioconda)
mamba create -c bioconda -c conda-forge -n pantas \
             python=3.10 biopython gffutils intervaltree bcftools samtools gffread vg=1.50.1 snakemake-minimal
mamba activate pantas
```

## pantas pipeline
``` sh
# Augment the annotated spliced pangenome with alignment information (run this for each replicate)
./pantas augment [condition1-rep1.gaf] [spliced-pangenome.annotated.gfa] > [condition1-rep1.gfa]

# Call events from each **augmented** graph
./pantas call [sample.gfa] [annotation.gtf] > [sample-events.csv]

# Quantify events across conditions (provide the two conditions with comma-separated path to the events csv)
./pantas quant condition1-rep1.csv,condition1-rep2.csv,condition1-rep3.csv \
             condition2-rep1.csv,condition2-rep2.csv,condition2-rep3.csv > [quantification.csv]
             
# Remap quantification on the linear reference using the annotation
./pantas remap -i [minimum intron size for a novel junction] [quantification.csv] [annotation.gtf] > [remap.csv]
```

### Event calling
The `call` mode of pantas provides several arguments that can be used to tweak the event calling:
```
  -w <INT>       Minimum read count (default: 3)
  -W <INT>       Minimum read count for annotated events (default: -1)
  -n             Call novel events (default: False)
  -a             Do not call known annotated events (default: False)
```

To call events from a reduced spliced pangenome, it is necessary to use the `-p` argument to provide the list of reference paths in the reduced graph (this file is created during the graph construction/indexing):
``` sh
python3 ./scripts/call.py --rp [spliced-pangenome.refpath] [sample.gfa] [annotation.gtf] > [sample-events.csv]
```

### Input preparation
The input of pantas are: an annotated spliced pangenome and the replicates aligned to this graph.

To build and index an annotated spliced pangenome, we provide a snakemake pipeline (`build/build.smk`): 
``` sh
snakemake -s build/build.smk -c4 --config fa=/path/to/reference.fa gtf=/path/to/annotation.gtf vcf=/path/to/variants.vcf.gz wd=/path/to/out/dir
```
Build step could also be called directlt from `pantas`:

``` sh
./pantas build -t [threads] -o [/path/to/out/dir] [reference.fa] [annotation.gtf] [variants.vcf.gz]
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

<!-- To build/index a **reduced** annotated spliced pangenomes, i.e., a graph representing a panel of genes of interest: -->
<!-- ``` sh -->
<!-- snakemake -s index-reduced.smk -c4 --config fa=/path/to/reference.fa gtf=/path/to/panel.gtf vcf=/path/to/variants.vcf.gz wd=/path/to/out/dir -->
<!-- ``` -->
<!-- The reduced annotated spliced pangenome and the index are stored in the `wd` directory: -->
<!-- ``` -->
<!-- # Annotated spliced pangenome in GFA format: -->
<!-- spliced-pangenes.annotated.gfa -->
<!-- # Compressed graph: -->
<!-- spliced-pangenes.xg -->
<!-- # Index: -->
<!-- spliced-pangenes.dist          -->
<!-- spliced-pangenes.gcsa -->
<!-- spliced-pangenes.gcsa.lcp -->
<!-- # Reduced reference paths: -->
<!-- spliced-pangenes.refpath -->
<!-- ``` -->

<!-- **Note:** using reduced annotated spliced pangenomes (when possible) is recommended since it hugely improves running times and RAM usage. -->

To map each replicate to the annotated spliced pangenome, we suggest to use `vg mpmap`:
``` sh
vg mpmap -x [spliced-pangenome.xg] -g [spliced-pangenome.gcsa] -d [spliced-pangenome.dist] -f [sample_1.fq] -f [sample_2.fq] -F GAF > [sample.gaf]
```

## Example
The `example` subdirectory contains example data that can be used to test pantas:
``` sh
# Prepare the graph
# This should take ~1 minute
snakemake -s build/build.smk -c4 --config fa=example/4.fa gtf=example/4.gtf vcf=example/4.vcf.gz wd=example/pantas-index

# Align the RNA-Seq sample to the graph
# This should take ~10 seconds
vg mpmap -x example/pantas-index/spliced-pangenome.xg \
         -g example/pantas-index/spliced-pangenome.gcsa \
	 -d example/pantas-index/spliced-pangenome.dist \
	 -f example/reads_1.fq -f example/reads_2.fq -F GAF > example/reads.gaf

# Augment the annotated spliced pangenome with alignment information
# This should be immediate
./pantas augment example/reads.gaf example/pantas-index/spliced-pangenome.annotated.gfa > example/reads.gfa

# Call all annotated events with minimum support 0 (since example RNA-Seq sample is very small)
# Note that using -W 0 is equivalent to extract all events from the graph
# This should take less than 2 seconds
./pantas call -W 0 example/reads.gfa example/4.gtf > example/reads.events.csv

# Quantify the events across the two conditions (an an example here we are using the same file twice)
# This should be immediate
./pantas quant example/reads.events.csv example/reads.events.csv > example/quant.csv

# Remap step
# This should be immediate
./pantas remap -i 25 example/quant.csv example/4.gtf > example/remap.csv
```

## Custom output format
#### Annotated and augmented spliced pangenome
The annotated spliced pangenome augmented with alignment information (output of `augment` mode of pantas) is stored in a GFA file where optional fields are used to store the annotation. We refer to the [documentation](docs/README.md).

#### Events
The events (output of `call` mode) are stored in a CSV file:
* event type (ES, A3, A5, IR)
* annotated/novel
* chromosome (e.g., 4)
* gene name (e.g., FBgn0004859)
* strand (e.g., +)
* junction1, based on annotation (e.g., FBtr0308074.4.5 or `?` if novel)
* junction1, in graph space (e.g., 2057>2065, meaning the junction link segments 2057 and 2065)
* junction1, on linear reference (e.g., 4:50614-50744)
* support for junction1 (e.g., 3)
* junction2, junction2 based on annotation, junction2 in graph space, junction2 on linear reference, support for junction2
* junction3, junction3 based on annotation, junction3 in graph space, junction3 on linear reference, support for junction3

We note that in the case of an exon skipping (or a cassette exon), all three junctions will be reported. In the case of an alternative splice site events, only two junctions are reported (1 and 2). A point (`.`) indicates that the junction is not used in the event.

#### Quantification
The differential quantification across conditions (output of `quant` mode of pantas) is stored in a CSV file:
* event type
* annotated/novel
* chromosome
* gene name
* strand
* junction1 (name on linear reference and node labels)
* junction2 (name on linear reference and node labels)
* junction3 (name on linear reference and node labels)
* support for canonical isoform involved in the event (one value per condition, separated by /)
* support for minor isoform involved in the event (one value per condition, separated by /)
* PSI value for condition 1
* PSI value for condition 2
* ΔPSI

#### Remapping
The differential quantification across conditions remapped on the linear reference (output of `remap` mode of pantas) is stored in a CSV file:
* event type
* annotated/novel
* refernce/haplotype
* chromosome
* gene name
* strand
* junction1 (name on linear reference, node labels, and position on linear reference)
* junction2 (name on linear reference, node labels, and position on linear reference)
* junction3 (name on linear reference, node labels, and position on linear reference)
* support for canonical isoform involved in the event (one value per condition, separated by /)
* support for minor isoform involved in the event (one value per condition, separated by /)
* PSI value for condition 1
* PSI value for condition 2
* ΔPSI


## Experiments
Experimental evaluation scripts can be found in the `./exps` subdirectory of this repository. We provide three snakemake pipelines which also contain more information on how to use pantas.
* `./exps/1-dm-sim/` is the evaluation on simulated data from Drosophila Melanogaster 
* `./exps/2-dm-real/` is the evaluation on real data from Drosophila Melanogaster
* `./exps/3-homo-real/` is the evaluation on real data from human

Additional details can be found in the README files available in these subdirectories.

## Authors
pantas is developed by [Simone Ciccolella](https://github.com/sciccolella), [Davide Cozzi](https://github.com/dlcgold), and [Luca Denti](https://github.com/ldenti).

For inquiries on this software please open an [issue](https://github.com/algolab/pantas/issues).
