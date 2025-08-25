#!/bin/bash

# default file extension .jpg
EXTENSION="jpg"

# declare empty array of positional arguments
declare -a POSITIONAL_ARGS

# define help string
HELP_STRING=$'Usage
    examples.bash [options] <source-file-names>
    echo | examples.bash [options] <source-file-names>

Specify the names of the source files to execute.

Options
    -h|--help                   Help message.
    -e|--extension  [jpg|eps]   File extension of output image. Defaults to jpg.\n'

# parse arugments
while [[ $# -gt 0 ]]; do
  case $1 in
        # file extension (.jpg or .eps)
    -e|--extension)
        EXTENSION="$2"
        shift
        shift
      ;;
    -h|--help)
        echo "$HELP_STRING"
        shift
        ;;
    -*|--*)
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

PATH_OUT=$PATH_WORKING_DIR/examples/cpp/out
mkdir -p $PATH_OUT

for FILE in "${POSITIONAL_ARGS[@]}"; do
    $PATH_EXAMPLE_BINS/$FILE $PATH_OUT/$FILE.$EXTENSION
done