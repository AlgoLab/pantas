```sh
git submodule update --init --recursive

cd deps/sdsl-lite
bash install.sh .

cd ../gbwt
bash install.sh .

cd ../libbdsg-easy
make -j4

cd ../..
g++ -O3 -o annotate annotate.cpp -I./deps/gbwt/include -I./deps/sdsl-lite/include -I./deps/libbdsg-easy/include/ -L./deps/sdsl-lite/lib -L./deps/gbwt/lib/ -L./deps/libbdsg-easy/lib/ -lgbwt -lsdsl -fopenmp -lbdsg -lhandlegraph

# Create the annotated spliced pangenome
snakemake -c [THREADS] --config fa=[FA] gtf=[GTF] vcf=[VCF] wd=[OUTDIR] -s index.smk
# graph is: [OUTDIR]/pantranscriptome.xg
# annotated graph is [OUTDIR]/pantranscriptome-annotated.gfa

# Index the graph
vg index --progress --threads [THREADS] --temp-dir [TMP-DIR] --gcsa-out pantranscriptome.gcsa --dist-name pantranscriptome.dist [OUTDIR]/pantranscriptome.xg
```