#!/bin/bash
mkdir -p benchmarks/logs
$BENCHMAKR_BINS/circleBenchmark --benchmark_out=$WORKING_DIRECTORY/benchmarks/logs/benchmark.json --benchmark_out_format=json --benchmark_counters_tabular=true --benchmark_time_unit=ms