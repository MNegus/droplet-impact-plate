#!/bin/bash

export OMP_NUM_THREADS=8

rm *.s
rm *.s.d
rm *.tests 
rm *.deps
rm *.tst

make uniform_plate.tst


