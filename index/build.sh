#!/bin/sh

set -e

bd=$(dirname $0)

FA=$1
GTF=$2
VCF=$3
WD=$4
NT=$5

samtools faidx $FA

mkdir -p $WD
mkdir -p $WD/times

for chrom in $(cut -f1 $FA.fai); do
	echo "--- $chrom ---"

	# This will be the inputs for the chromosome
	fa=$WD/$chrom.fa
	gtf=$WD/$chrom.gtf
	vcf=$WD/$chrom.vcf.gz
	wd=$WD/$chrom

	# Prepare input
	samtools faidx $FA $chrom >$fa
	samtools faidx $fa
	grep -P "^$chrom\t" $GTF >$gtf
	if [ ! -f $vcf ]; then
		bcftools view -Oz $VCF $chrom >$vcf
		tabix -p vcf $vcf
	fi

	mkdir -p $wd
	td=$WD/times/$chrom
	mkdir -p $td

	if [ ! -f $wd/pangenome_walts.pg ]; then
		echo "[ $(date) ] $chrom - vg construct with alts"
		/usr/bin/time -vo $td/1-construct.time sh -c "vg construct --threads $NT -r $fa -v $vcf --alt-paths --node-max 1024 | vg convert --packed-out - > $wd/pangenome_walts.pg"
	fi
	if [ ! -f $wd/spliced_pangenome_walts.pg ]; then
		echo "[ $(date) ] $chrom - vg rna (first pass)"
		/usr/bin/time -vo $td/2-rna.time vg rna --threads $NT --transcripts $gtf $wd/pangenome_walts.pg >$wd/spliced_pangenome_walts.pg
	fi

	if [ ! -f $wd/spliced_pangenome.pg ]; then
		echo "[ $(date) ] $chrom - vg drop alts"
		/usr/bin/time -vo $td/3-drop.time sh -c "vg paths --drop-paths --variant-paths -x $wd/spliced_pangenome_walts.pg | vg convert --packed-out - > $wd/spliced_pangenome.pg"
	fi

	if [ ! -f $wd/reference.gbwt ]; then
		echo "[ $(date) ] $chrom - vg gbwt (reference path)"
		/usr/bin/time -vo $td/4-gbwt.time vg gbwt --index-paths --num-jobs $NT --xg-name $wd/spliced_pangenome.pg --output $wd/reference.gbwt
	fi

	if [[ $(bcftools view -H -G $vcf | head) ]]; then
		if [ ! -f $wd/haplotypes.gbwt ]; then
			echo "[ $(date) ] $chrom - vg gbwt (haplotypes)"
			/usr/bin/time -vo $td/5a-gbwt.time vg gbwt --discard-overlaps --num-jobs $NT --preset 1000gp --vcf-input $vcf --xg-name $wd/spliced_pangenome_walts.pg --output $wd/haplotypes.gbwt
		fi
		if [ ! -f $wd/samples.gbwt ]; then
			echo "[ $(date) ] $chrom - vg gbwt (merge)"
			/usr/bin/time -vo $td/5b-gbwt.time vg gbwt --merge -o $wd/samples.gbwt $wd/haplotypes.gbwt $wd/reference.gbwt
		fi
	else
		if [ ! -f $wd/samples.gbwt ]; then
			echo "[ $(date) ] $chrom - no haplotypes. simply copying"
			/usr/bin/time -vo $td/5a-gbwt.time cp $wd/reference.gbwt $wd/samples.gbwt
		fi
	fi

	if [ ! -f $wd/pantranscriptome.pg ]; then
		echo "[ $(date) ] $chrom - vg rna (second pass) - log: $wd/pantranscriptome.log"
		# --add-ref-paths and --add-hap-paths are needed since we need to prune but keep transcripts
		/usr/bin/time -vo $td/6-rna.time vg rna --progress --threads $NT --add-ref-paths --add-hap-paths --haplotypes $wd/samples.gbwt --transcripts $gtf --write-info $wd/pantranscriptome.info --write-gbwt $wd/pantranscriptome.gbwt $wd/spliced_pangenome.pg >$wd/pantranscriptome.pg 2>$wd/pantranscriptome.log
	fi

	if [ ! -f $wd/pantranscriptome-pruned.pg ]; then
		echo "[ $(date) ] $chrom - vg prune"
		/usr/bin/time -vo $td/7-prune.time vg prune --restore-paths --progress --threads $NT $wd/pantranscriptome.pg >$wd/pantranscriptome-pruned.pg
	fi

	if [ ! -f $wd/pantranscriptome-pruned.gfa ]; then
		echo "[ $(date) ] $chrom - vg convert to gfa"
		/usr/bin/time -vo $td/8-view.time vg view --threads $NT $wd/pantranscriptome-pruned.pg >$wd/pantranscriptome-pruned.gfa
	fi

	if [ ! -f $wd/pantranscriptome-annotated.gfa ]; then
		echo "[ $(date) ] $chrom - annotate - log: $wd/annotation.log"
		/usr/bin/time -vo $td/9-annotate.time $bd/annotate $wd/pantranscriptome-pruned.gfa $wd/pantranscriptome.info $wd/samples.gbwt $wd/pantranscriptome.gbwt >$wd/pantranscriptome-annotated.gfa 2>$wd/annotation.log
	fi

	echo "[ $(date) ] $chrom - Done"
done

echo "[ $(date) ] combine"
/usr/bin/time -vo $WD/times/combine.time python3 $bd/combine.py $WD/*/pantranscriptome-annotated.gfa >$WD/pantranscriptome-annotated.gfa

echo "[ $(date) ] vg convert"
/usr/bin/time -vo $WD/times/convert.time vg convert --gfa-in --xg-out $WD/pantranscriptome-annotated.gfa >$WD/pantranscriptome.xg

# # We do not need paths in the graph. Tried to index with and without paths. Same gcsa and lcp (diff) and dist (vg view --distance-in)
# echo "[ $(date) ] vg index - log: $WD/index.log"
# /usr/bin/time -vo $WD/times/index.time vg index --progress --threads $NT --temp-dir $WD --gcsa-out $WD/pantranscriptome.gcsa --dist-name $WD/pantranscriptome.dist $WD/pantranscriptome.xg 2>$WD/index.log

echo "[ $(date) ] Done"
