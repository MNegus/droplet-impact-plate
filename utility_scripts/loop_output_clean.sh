#!/bin/bash

parent_dir=$1
code_name=$2

for MAXLEVEL in 10 11 12 13 14
do
    echo Max level = $MAXLEVEL
    ./output_clean.sh $parent_dir/max_level_$MAXLEVEL
done

