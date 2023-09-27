#!/bin/sh

mode=$1
ROOT=$(dirname $0)
echo $ROOT

if [[ $mode == "index" ]]
then
    FA=$2
    GTF=$3
    VCF=$4
    WD=$5
    NT=$6
    pushd $ROOT
    snakemake -s index.smk -c$NT --config fa=$FA gtf=$GTF vcf=$VCF wd=$WD
    popd
elif [[ $mode == "align" ]]
then
    INDEX=$2
    FQLIST=$3
    WD=$4
    NT=$5
    snakemake -s align.smk -c$NT --config index=$INDEX fqlist=$FQLIST wd=$WD
elif [[ $mode == "quantify" ]]
then
    echo "Quantify mode."
else
    echo "Unknown mode."
fi
