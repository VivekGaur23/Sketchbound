#!/bin/sh
set -e
CXX=${CXX:-g++}
PYTHON=${PYTHON:-python}
FLAGS="-std=c++20 -O2 -Wall -Iinclude -Ibench"
mkdir -p results figures
$CXX $FLAGS -o tests/test_sketches tests/test_sketches.cpp
$CXX $FLAGS -o bench/experiments bench/experiments.cpp
./tests/test_sketches
./bench/experiments
$PYTHON analysis/plot.py
