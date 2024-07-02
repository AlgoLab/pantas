#!/bin/sh

fa=$1
gtf=$2

grep -P "\tgene\t" $gtf | while read line
do
    c=$(echo $line | cut -f1 -d' ')
    s=$(echo $line | cut -f4 -d' ')
    e=$(echo $line | cut -f5 -d' ')
    gidx=$(echo $line | grep -E -o 'gene_id "[0-9A-Za-z]+";' | cut -f2 -d'"')
    echo ">$gidx $c:$s-$e"
    samtools faidx $fa $c:$s-$e | sed '1d'
done
