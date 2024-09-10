# Exps on simulated data from Drosophila

### Setup input
``` sh
wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/fasta/dmel-all-chromosome-r6.51.fasta.gz
gunzip dmel-all-chromosome-r6.51.fasta.gz
samtools faidx dmel-all-chromosome-r6.51.fasta 2L 2R 3L 3R 4 X Y > ref.chroms.fa

wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/gtf/dmel-all-r6.51.gtf.gz
gunzip dmel-all-r6.51.gtf.gz
# gene FBgn0002781 has inverted transcripts
for c in 2L 2R 3L 3R 4 X Y ; do grep -P "^$c\t" dmel-all-r6.51.gtf ; done | grep -v FBgn0002781 > genes.chroms.gtf
python3 ./scripts/clean_gtf.py genes.chroms.gtf > genes.gtf
grep -P "\tgene\t" genes.gtf | cut -f1,4,5 > genes.bed

wget https://resources.aertslab.org/DGRP2/BCM-HGSC/final/dm6/DGRP2.source_BCM-HGSC.dm6.final.SNPs_only.vcf.gz.tbi
tabix -p vcf DGRP2.source_BCM-HGSC.dm6.final.SNPs_only.vcf.gz
bcftools view -R genes.chroms.bed -c 1 -q 0.01 -Oz DGRP2.source_BCM-HGSC.dm6.final.SNPs_only.vcf.gz > DGRP2.SNPs_only.genes.vcf.gz
tabix -p vcf DGRP2.SNPs_only.genes.vcf.gz
python3 ./scripts/fix_vidx.py DGRP2.SNPs_only.genes.vcf.gz | bcftools norm -d all -Oz -f ref.chroms.fa > DGRP2.SNPs_only.genes.norm.vcf.gz
tabix -p vcf DGRP2.SNPs_only.genes.norm.vcf.gz

wget https://resources.aertslab.org/DGRP2/BCM-HGSC/final/dm6/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz
tabix -p vcf DGRP2.source_BCM-HGSC.dm6.final.vcf.gz
bcftools view -R genes.bed -c 1 -q 0.01 -Oz DGRP2.source_BCM-HGSC.dm6.final.vcf.gz > DGRP2.genes.vcf.gz
tabix -p vcf DGRP2.genes.vcf.gz
python3 ./scripts/fix_vidx.py DGRP2.genes.vcf.gz | bcftools norm -d all -Oz -f ref.chroms.fa > DGRP2.genes.norm.vcf.gz
tabix -p vcf DGRP2.genes.norm.vcf.gz
```


## Setup environment
``` sh
mamba create -n pantas-simexps -c bioconda -c conda-forge python=3.11.5 r-base numpy samtools bcftools biopython intervaltree snakemake-minimal vg=1.56
conda activate pantas-simexps
R
> install.packages("remotes")
> remotes::install_github("biomedbigdata/ASimulatoR")
> q()
```

### Run experiments
```
snakemake -pc16 --config fa=/path/to/reference.fa gtf=/path/to/annotation.gtf vcf=/path/to/population.vcf.gz odir=/path/to/output/directory n=<nreads> --use-conda [-n]
# Results: $odir/*/compare-{anno,novel}.csv
python3 ./scripts/plot_pr.py $odir/*/compare-anno.csv
python3 ./scripts/plot_pr.py $odir/*/compare-novel.csv
```
