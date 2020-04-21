#!/bin/bash

# Deletes the previous directory and test files
rm -r ode


# Runs the code using the Makefile
for RUN_NO in 0 1 2 3 4 5 6 7 8 9 10 11 12 13
do
    rm *.s
    rm *.s.d
    rm *.tests 
    rm *.deps
    rm *.tst

    make ode.tst RUN_NO=$RUN_NO 
done

# Moves the data
mv -r ode ../raw_data
