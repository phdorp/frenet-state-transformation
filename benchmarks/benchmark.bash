#!/bin/bash

# default file extension
EXTENSION="json"

# default time unit
TIME_UNIT="ms"

# default tabular counters
COUNTERS_TABULAR="true"

# declare empty array of positional arguments
declare -a POSITIONAL_ARGS

# declare empty array of additional arguments
declare -a ADDITIONAL_ARGS

# current datetime
DATETIME=$(eval date +"%y%m%d-%H%M%S")

# define help string
HELP_STRING=$'Usage
    benchmark.bash [options] <source-file-names>

Specify the names of the source files to execute.
Additional options in the format --option=value are forwared to the executable.
Look at https://github.com/google/benchmark/blob/main/docs/user_guide.md#running-a-subset-of-benchmarks for more details.

Options
    -h|--help                                           Help message.
    --benchmark_out_format         [json|console|csv]   Write result to file or console. Defaults to json.
    --benchmark_time_unit          [ns|us|ms|s]         Time unit for time measurements. Defaults to ms.
    --benchmark_counters_tabular   [true|false]         Whether to print costum counters in tabular format. Defaults to true.\n'

# parse arugments
while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
        echo "$HELP_STRING"
        shift
        ;;
    --benchmark_out_format)
        EXTENSION="$2"
        shift
        shift
        ;;
    --benchmark_time_unit)
        TIME_UNIT="$2"
        shift
        shift
        ;;
    --benchmark_counters_tabular)
        COUNTERS_TABULAR="$2"
        shift
        shift
        ;;
    --*)
        ADDITIONAL_ARGS+=("$1")
        shift
        shift
        ;;
    -*)
        echo "Unknown option $1"
        exit 1
        ;;
    # executables as positional arguments
    *)
        POSITIONAL_ARGS+=("$1")
        shift
        ;;
  esac
done

OUT_DIRECTORY=$WORKING_DIRECTORY/benchmarks/logs
mkdir -p $OUT_DIRECTORY

for FILE in "${POSITIONAL_ARGS[@]}"; do
    $BENCHMARK_BINS/$FILE --benchmark_out=$OUT_DIRECTORY/"$DATETIME"_"$FILE".$EXTENSION --benchmark_out_format=$EXTENSION --benchmark_counters_tabular=$COUNTERS_TABULAR --benchmark_time_unit=$TIME_UNIT "${ADDITIONAL_ARGS[@]}"
done