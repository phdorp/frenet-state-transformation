#!/bin/bash
for example in $*;
do
    OUT_DIRECTORY=$WORKING_DIRECTORY/examples/out
    mkdir -p OUT_DIRECTORY
    $EXAMPLE_BINS/$example $OUT_DIRECTORY/$example.jpg
done