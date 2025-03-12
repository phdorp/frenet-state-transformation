#!/bin/bash
mkdir -p benchmarks/logs
for commit in $*;
do
    git checkout $commit
    rm -r build <&-
    mkdir build
    cmake -B build -D BUILD_BENCHMARK=ON -D CMAKE_CXX_COMPILER=g++-13
    cmake --build build
    DATETIME=$(git show -s --format=%cd --date=format:'%y%m%d-%H%M' $commit)
    $BENCHMAKR_BINS/circleBenchmark --benchmark_out=$WORKING_DIRECTORY/benchmarks/logs/$DATETIME-$commit.json --benchmark_out_format=json --benchmark_counters_tabular=true --benchmark_time_unit=ms --benchmark_context=commit=$commit
done