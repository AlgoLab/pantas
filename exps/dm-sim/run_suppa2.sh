#!/usr/bin/env sh

#!/usr/bin/env sh

set -xe

ioe=$1
inquants=$2
outd=$3

wd=$(dirname $0)

mkdir -p $outd
multipleFieldSelection.py -i $inquants -k 1 -f 4 -o $outd/all.iso_tpm.txt
suppa.py psiPerEvent -i $ioe -e $outd/all.iso_tpm.txt -o $outd/events
cols=$(head -n1 $outd/events.psi)
IFS=$'\t' read -ra h_tot <<< "$cols"
h1=${h_tot[@]:0:$((${#h_tot[@]}/2))}
h2=${h_tot[@]:$((${#h_tot[@]}/2)):${#h_tot[@]}/2}
h1s=$(echo ${h1[@]} | tr ' ' ',')
h2s=$(echo ${h2[@]} | tr ' ' ',')
echo $ioe
echo $inquants
echo $outd
echo $wd
echo $h1s
echo $h2s
Rscript $wd/split_file.R $outd/all.iso_tpm.txt $h1s $h2s $outd/h1.tpm $outd/h2.tpm
Rscript $wd/split_file.R $outd/events.psi $h1s $h2s $outd/h1.psi $outd/h2.psi
suppa.py diffSplice -m empirical -i $ioe -e $outd/h1.tpm $outd/h2.tpm -p $outd/h1.psi $outd/h2.psi -o $outd/DIFF
if [ -f $outd/DIFF.dpsi.temp.0 ]; then
    mv $outd/DIFF.dpsi.temp.0 $outd/DIFF.dpsi
fi
