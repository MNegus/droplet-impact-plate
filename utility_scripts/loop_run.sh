#!/bin/bash

parent_dir=$1
code_name=$2
move_dest=$3

for MAXLEVEL in 10 11 12 13 14
do
    echo Max level = $MAXLEVEL
    cd $parent_dir/max_level_$MAXLEVEL/code

    ./run_simulation.sh $code_name 16
done

