#!/bin/bash

# ./output_clean.sh /mnt/newarre/cantilever_paper_data/alpha_varying/alpha_2
# ./output_clean.sh /mnt/newarre/cantilever_paper_data/alpha_varying/alpha_5
# ./output_clean.sh /mnt/newarre/cantilever_paper_data/alpha_varying/alpha_10
# ./output_clean.sh /mnt/newarre/cantilever_paper_data/alpha_varying/alpha_20
# ./output_clean.sh /mnt/newarre/cantilever_paper_data/alpha_varying/alpha_100


for LEVEL in 10 11 12 13 14
do  
    ./output_clean.sh /media/michael/newarre/supplementary_material/validation/stationary_validation/level_$LEVEL
    ./output_clean.sh /media/michael/newarre/supplementary_material/validation/moving_validation/level_$LEVEL
done
