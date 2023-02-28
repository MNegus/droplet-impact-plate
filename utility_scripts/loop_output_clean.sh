#!/bin/bash

parent_dir=$1

# for MAXLEVEL in 10 11 12 13 14
# do
#     echo Max level = $MAXLEVEL
#     ./output_clean.sh $parent_dir/max_level_$MAXLEVEL
# done
for DIR_NAME in $parent_dir/*/
do
    echo $DIR_NAME
    ./output_clean.sh $DIR_NAME
done

