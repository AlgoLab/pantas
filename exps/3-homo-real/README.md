# Exps on real data from Human

### Setup input
``` sh
# Reference genome
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa $(seq 1 22) X Y > reference.chroms.fa
sed -i "s/^>/>chr/g" reference.chroms.fa
samtools faidx reference.chroms.fa

# Gene Annotation
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
for c in $(seq 1 22) X Y ; do zgrep -P "^$c\t" Homo_sapiens.GRCh38.109.gtf.gz ; done | sed -e 's/^/chr/' > annotation.gtf

# Variants
mkdir 1kgp-GRCh38
cd 1kgp-GRCh38
# get all VCF from ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
for c in $(seq 1 22) X ; do echo $c ; bcftools view -c 1 -q 0.001 1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -Oz | bcftools norm -d all -Oz -f ../reference.chroms.fa > chr$c-mod.vcf.gz ; tabix -p vcf chr$c-mod.vcf.gz ; done
bcftools concat -Oz chr*-mod.vcf.gz > ../1kgp-GRCh38.vcf.gz
tabix -p vcf ../1kgp-GRCh38.vcf.gz

# RNA-Seq data
# fasterq-dump SRX651011
for fq in $(ls *.fastq) ; do bn=$(basename $fq .fastq) ; echo $bn ; /usr/bin/python3 clean_reads_from_N.py $bn.fastq > $bn.clean.fq ; done

# Get truth
# truth-pos.csv is Table S4 from SUPPA2 paper
# truth-neg.csv is Table S10 from SUPPA2 paper
cd data
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
cut -f1 -d';' truth-pos.csv | cut -d',' -f 2- | awk -F, '{print $2,$3,$4,$7,$6}' OFS=,  | tr ',' '\t' > truth-pos.bed
./liftOver truth-pos.bed hg19ToHg38.over.chain truth-pos.hg38.bed truth-pos.unlifted.bed
cut -d',' -f2- truth-neg.csv | tr ',' '\t' | tail -n +2 > truth-neg.bed
./liftOver truth-neg.bed hg19ToHg38.over.chain truth-neg.hg38.bed truth-neg.unlifted.bed
cut -f4 truth-pos.hg38.bed truth-neg.hg38.bed | sort -u | while read idx ; do echo \"$idx\"; done > genes.list
# Some genes cannot be found in the GTF due to different name, + 3 true events have wrong coordinates after liftover
while read old new ; do sed -i "s/$old/$new/g" genes.list ; done < gene-remapping.txt
cat <(sed -e "s/^/POS\t/" truth-pos.hg38.bed) <(sed -e "s/^/NEG\t/" truth-neg.hg38.bed) | grep -P -v "MPHOSPH10|ZCCH9C|ZZZ3" > truth.tsv
```

### Run experiments
``` sh
# Update config/config.yaml
snakemake -pj 32 --use-conda [-n]

# Analysis and plot, assuming WD to be the odir
python3 ./workflow/scripts/compare.py truth.tsv $WD/pantas2/quant.w5.csv $WD/rMATS/SE.MATS.JC.txt $WD/whippet/psi.diff $WD/suppa2/OUT/DIFF.dpsi

# Check if true events are annotated/novel
python3 ./workflow/scripts/check_novel.py data/truth.tsv genes.gtf

# Check if events are covered by alignments in the provided BAMs (missed.txt is part of compare.py output)
python3 ./workflow/scripts/check_coverage.py data/missed.txt genes.gtf SRR1513329.bam SRR1513330.bam SRR1513331.bam SRR1513332.bam SRR1513333.bam SRR1513334.bam
```

