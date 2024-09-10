#!/bin/sh

set -e

git submodule update --init --recursive

cd build/deps/sdsl-lite
bash install.sh .

cd ../gbwt
bash install.sh .

cd ../libbdsg-easy
make -j4

cd ../..
g++ -O3 -o annotate annotate.cpp -I./deps/gbwt/include -I./deps/sdsl-lite/include -I./deps/libbdsg-easy/include/ -L./deps/sdsl-lite/lib -L./deps/gbwt/lib/ -L./deps/libbdsg-easy/lib/ -lgbwt -lsdsl -fopenmp -lbdsg -lhandlegraph

cd ..

echo "--- Everything done! ---"
