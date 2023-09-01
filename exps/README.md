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
```
