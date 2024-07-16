```sh
git submodule update --init --recursive

cd deps/sdsl-lite
bash install.sh .

cd ../gbwt
bash install.sh .

# this is not needed anymore
# cd ../libbdsg-easy
# make -j4

cd ../..
g++ -O3 -o annotate annotate.cpp -I./deps/gbwt/include -I./deps/sdsl-lite/include -L./deps/sdsl-lite/lib -L./deps/gbwt/lib/ -lgbwt -lsdsl -fopenmp

# Create the annotated spliced pangenome
bash build.sh [FA] [GTF] [VCF] [OUTDIR] [THREADS]
# graph is: [OUTDIR]/pantranscriptome.xg
# annotated graph is [OUTDIR]/pantranscriptome-annotated.gfa

# Index the graph
# Note: we do not need paths in the graph. I've tried to index with and without paths.
# Same .gcsa and .lcp (diff) and .dist (vg view --distance-in)
vg index --progress --threads [THREADS] --temp-dir [TMP-DIR] --gcsa-out pantranscriptome.gcsa --dist-name pantranscriptome.dist [OUTDIR]/pantranscriptome.xg
```