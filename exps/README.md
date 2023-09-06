### Drosophila

``` sh
mkdir chroms

wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/fasta/dmel-all-chromosome-r6.51.fasta.gz
gunzip dmel-all-chromosome-r6.51.fasta.gz
samtools faidx dmel-all-chromosome-r6.51.fasta 2L 2R 3L 3R 4 X Y > ref.chroms.fa


wget http://ftp.flybase.net/releases/FB2023_02/dmel_r6.51/gtf/dmel-all-r6.51.gtf.gz
gunzip dmel-all-r6.51.gtf.gz
for c in 2L 2R 3L 3R 4 X Y ; do grep -P "^$c\t" dmel-all-r6.51.gtf ; done > genes.chroms.gtf
grep -P "\tgene\t" genes.chroms.gtf | cut -f1,4,5 > genes.chroms.bed

wget https://resources.aertslab.org/DGRP2/BCM-HGSC/final/dm6/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz
tabix -p vcf DGRP2.source_BCM-HGSC.dm6.final.vcf.gz
bcftools view -R genes.chroms.bed -c 1 -q 0.01 -Oz DGRP2.source_BCM-HGSC.dm6.final.vcf.gz > DGRP2.genes.vcf.gz
tabix -p vcf DGRP2.genes.vcf.gz
python3 fix_vidx.py DGRP2.genes.vcf.gz | bcftools norm -d all -Oz -f ref.chroms.fa > DGRP2.genes.norm.vcf.gz
tabix -p vcf DGRP2.genes.norm.vcf.gz

# check GTF file
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir dm6.genome-STAR --genomeFastaFiles ref.chroms.fa --sjdbGTFfile genes.chroms.mod.gtf --sjdbOverhang 99

```

#### Simulate
``` sh
Rscript exps/asimulator.R INPUT OUTPUT READ_LEN SEQ_DEPTH
# eg: Rscript exps/asimulator.R /data/drosmel-age-recomb/simulated/input /data/drosmel-age-recomb/simulated/output-3 150 2000000

python exps/filter_reads.py OUTPU/sample_N_1.fastq OUTPUT/sample_N_2.fastq
# eg: python exps/filter_reads.py /data/drosmel-age-recomb/simulated/output-3/sample_01_1.fastq /data/drosmel-age-recomb/simulated/output-3/sample_01_2.fastq

# No need to use sample_N_2 the headers are the same
python exps/simrc.py OUTPUT/sample_N_1.clean.fq OUTPUT/exon_junction_coverage.tsv > OUTPUT/readcount.sample_N.csv
# eg: python exps/simrc.py /data/drosmel-age-recomb/simulated/output-3/sample_01_1.clean.fq /data/drosmel-age-recomb/simulated/output-3/exon_junction_coverage.tsv

STAR --runThreadN 16 --genomeDir dm6.genome-STAR \
        --readFilesIn OUTPUT/sample_N_1.clean.fq OUTPUT/sample_N_2.clean.fq \
        --outFileNamePrefix OUTPUT/STAR/sample_N \
        --outSAMtype BAM SortedByCoordinate --outSAMattributes All

# eg: STAR --runThreadN 16 --genomeDir /data/drosmel-age-recomb/dm6.genome-STAR \
#         --readFilesIn /data/drosmel-age-recomb/simulated/output-3/sample_01_1.clean.fq /data/drosmel-age-recomb/simulated/output-3/sample_01_2.clean.fq \
#         --outFileNamePrefix /data/drosmel-age-recomb/simulated/output-3/STAR/sample_01 \
#         --outSAMtype BAM SortedByCoordinate --outSAMattributes All

```