#!/bin/sh

inputs=( "$@" )
echo "Tool,"$(head -1 $inputs)

for f in "${@}"
do
    bn=$(basename $f .res.csv)
    tail -n+2 $f | while read line
    do
	echo $bn,$line
    done
done
