# pantas

This repository contains the codebase of pantas, a pangenomic approach for performing haplotype-aware differential AS events quantification across RNA-Seq conditions. pantas is based on the notion of annotated spliced pangenomes, which are [spliced pangenomes](https://doi.org/10.1038/s41592-022-01731-9) augmented with additional information needed for AS events inference and quantification.

Alongside pantas, we provide a set of utilities to build and index a (annotated) spliced pangenome. We devised [two alternative construction pipelines](https://github.com/AlgoLab/pantas/tree/main#input-preparation), one for full genome analysis and one for the analysis of a panel of genes of interest (reduced indexing).


## Installation
``` sh
git clone https://github.com/AlgoLab/pantas.git

# Install dependencies (all available from bioconda)
mamba create -c bioconda -c conda-forge -n pantas \
             python=3.10 biopython gffutils intervaltree bcftools samtools gffread vg=1.56 snakemake-minimal
mamba activate pantas

cd pantas
bash setup.sh
# if setup ends properly, you should see the message "--- Everything done! ---"
```

## pantas pipeline
All steps of pantas can be easily run using the `pantas` script.

``` sh
# Print running modes
./pantas -h
```

All modes print a help message if run with the `-h` argument. Here we provide the minimum command lines arguments that are required to run `pantas`:
``` sh
# Augment the annotated spliced pangenome with alignment information (run this for each replicate)
./pantas augment [condition1-rep1.gaf] [spliced-pangenome.annotated.gfa] > [condition1-rep1.gfa]

# Call events from each **augmented** graph
./pantas call [sample.gfa] [annotation.gtf] > [sample-events.csv]

# Quantify events across conditions (provide the two conditions with comma-separated path to the events csv)
./pantas quant condition1-rep1.csv,condition1-rep2.csv,condition1-rep3.csv \
             condition2-rep1.csv,condition2-rep2.csv,condition2-rep3.csv > [quantification.csv]
             
# Remap quantification on the linear reference using the annotation
./pantas remap [quantification.csv] [annotation.gtf] > [quantification.remap.csv]
```

### Input preparation
The input of pantas are: an annotated spliced pangenome and the replicates aligned to this graph.

To build and index an annotated spliced pangenome, we provide a snakemake pipeline (`build/build.smk`) that can be conveniently run via the `pantas` script: 

``` sh
./pantas build -t [threads] -o [/path/to/out/dir] [reference.fa] [annotation.gtf] [variants.vcf.gz]
```

The annotated spliced pangenome is stored in the working directory directory (`-o` argument):
``` sh
# Annotated spliced pangenome in GFA format:
pantranscriptome-annotated.gfa
# Compressed graph:
pantranscriptome.xg

The graph can then be indexed using `vg index`:
``` sh
vg index --gcsa-out [pantranscriptome.gcsa] --dist-name [pantranscriptome.dist] [pantranscriptome.xg]
```

This will produce 3 files in the working directory:
``` sh
# Index:
pantranscriptome.dist
pantranscriptome.gcsa
pantranscriptome.gcsa.lcp
```

To map each replicate to the annotated spliced pangenome, we suggest to use `vg mpmap`:
``` sh
vg mpmap -x [pantranscriptome.xg] -g [pantranscriptome.gcsa] -d [pantranscriptome.dist] -f [sample_1.fq] -f [sample_2.fq] -F GAF > [sample.gaf]
```

## Example
The `example` subdirectory contains example data that can be used to test pantas:
``` sh
# Prepare the graph
./pantas build -t 4 -o example/pantas-index example/4.fa example/4.gtf example/4.vcf.gz

# Index the graph
vg index --progress --threads 4 --gcsa-out example/pantas-index/pantranscriptome.gcsa --dist-name example/pantas-index/pantranscriptome.dist example/pantas-index/pantranscriptome.xg

# Align the RNA-Seq sample to the graph
vg mpmap -x example/pantas-index/pantranscriptome.xg \
         -g example/pantas-index/pantranscriptome.gcsa \
         -d example/pantas-index/pantranscriptome.dist \
         -f example/reads_1.fq -f example/reads_2.fq -F GAF > example/reads.gaf

# Augment the annotated spliced pangenome with alignment information
./pantas augment example/reads.gaf example/pantas-index/pantranscriptome-annotated.gfa > example/pantranscriptome-annotated-wreads.gfa

# Call all annotated events with minimum support 0 (since example RNA-Seq sample is very small)
# Note that using -w 0 is equivalent to extract all events from the graph
./pantas call -w 0 example/pantranscriptome-annotated-wreads.gfa example/4.gtf > example/reads.events.csv

# Quantify the events across the two conditions (an an example here we are using the same file twice)
./pantas quant example/reads.events.csv example/reads.events.csv > example/quant.csv

# Remap step
./pantas remap example/quant.csv example/4.gtf > example/quant-remap.csv
# this should produce 205 events
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
* transcripts annotation for junctions 1, 2, and 3 (3 fields)
* edges for junctions 1, 2, and 3 (3 fields)
* support for canonical isoform involved in the event (one value per condition, separated by /)
* support for minor isoform involved in the event (one value per condition, separated by /)
* PSI value for condition 1
* PSI value for condition 2
* ΔPSI

#### Remapping
The differential quantification across conditions remapped on the linear reference (output of `remap` mode of pantas) is stored in a CSV file:
* event type
* annotated/novel
* reference/haplotype (telling if event involves reference transcripts or haplotype transcripts)
* chromosome
* gene name
* strand
* transcripts annotation for junctions 1, 2, and 3 (3 fields)
* edges for junctions 1, 2, and 3 (3 fields)
* reference coordinates for junctions 1, 2, and 3 (3 fields)
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
