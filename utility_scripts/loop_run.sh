#!/bin/bash
# Sequentially runs all of the scripts in a given parent directory

parent_dir=$1
code_name=$2
cores=$3

for DIR_NAME in $parent_dir/*/
do
    echo $DIR_NAME
    cd $DIR_NAME/code

    ./run_simulation.sh $code_name $cores
done
