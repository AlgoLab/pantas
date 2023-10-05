# Experiments
## Simulated on Drosophila

``` sh
wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/fasta/dmel-all-chromosome-r6.51.fasta.gz
gunzip dmel-all-chromosome-r6.51.fasta.gz
samtools faidx dmel-all-chromosome-r6.51.fasta 2L 2R 3L 3R 4 X Y > ref.chroms.fa


wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/gtf/dmel-all-r6.51.gtf.gz
gunzip dmel-all-r6.51.gtf.gz
# gene FBgn0002781 has inverted transcripts
for c in 2L 2R 3L 3R 4 X Y ; do grep -P "^$c\t" dmel-all-r6.51.gtf ; done | grep -v FBgn0002781 > genes.chroms.gtf
python3 ~/code/pantas2/exps/dm-sim/clean_gtf.py genes.chroms.gtf > genes.gtf
grep -P "\tgene\t" genes.gtf | cut -f1,4,5 > genes.bed


wget https://resources.aertslab.org/DGRP2/BCM-HGSC/final/dm6/DGRP2.source_BCM-HGSC.dm6.final.SNPs_only.vcf.gz.tbi
tabix -p vcf DGRP2.source_BCM-HGSC.dm6.final.SNPs_only.vcf.gz
bcftools view -R genes.chroms.bed -c 1 -q 0.01 -Oz DGRP2.source_BCM-HGSC.dm6.final.SNPs_only.vcf.gz > DGRP2.SNPs_only.genes.vcf.gz
tabix -p vcf DGRP2.SNPs_only.genes.vcf.gz
python3 fix_vidx.py DGRP2.SNPs_only.genes.vcf.gz | bcftools norm -d all -Oz -f ref.chroms.fa > DGRP2.SNPs_only.genes.norm.vcf.gz
tabix -p vcf DGRP2.SNPs_only.genes.norm.vcf.gz

wget https://resources.aertslab.org/DGRP2/BCM-HGSC/final/dm6/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz
tabix -p vcf DGRP2.source_BCM-HGSC.dm6.final.vcf.gz
bcftools view -R genes.chroms.bed -c 1 -q 0.01 -Oz DGRP2.source_BCM-HGSC.dm6.final.vcf.gz > DGRP2.genes.vcf.gz
tabix -p vcf DGRP2.genes.vcf.gz
python3 fix_vidx.py DGRP2.genes.vcf.gz | bcftools norm -d all -Oz -f ref.chroms.fa > DGRP2.genes.norm.vcf.gz
tabix -p vcf DGRP2.genes.norm.vcf.gz
```


## Simulation
``` sh
mamba create -n pantas2 -c bioconda -c conda-forge r-base samtools bcftools biopython intervaltree
conda activate pantas2
R
> install.packages("remotes")
> remotes::install_github("biomedbigdata/ASimulatoR")
> q()

cd dm-sim
snakemake -c16 --config fa=/path/to/reference.fa gtf=/path/to/annotation.gtf vcf=/path/to/population.vcf.gz odir=/path/to/output/directory n=<nreads>
```

## Evaluation
``` sh

```



### Human

``` sh
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa $(seq 1 22) X Y > reference.chroms.fa

wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
for c in $(seq 1 22) X Y ; do zgrep -P "^$c\t" Homo_sapiens.GRCh38.109.gtf.gz ; done > annotation.chroms.gtf

mkdir 1kgp-GRCh38
cd 1kgp-GRCh38
# get all VCF from ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

for chr in $(seq 1 22) X ; do bcftools view -c 1 -q 0.001 1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -Oz | python3 ~/code/pantas2/exps/fix_vidx.py | bcftools norm -d all -Oz -f ../reference.chroms.fa > $c-mod.vcf.gz ; tabix -p vcf $c-mod.vcf.gz ; done
bcftools concat -Oz *-mod.vcf.gz > ../1kgp-GRCh38.vcf.gz
tabix -p vcf ../1kgp-GRCh38.vcf.gz
```
