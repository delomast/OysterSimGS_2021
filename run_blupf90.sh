#!/bin/bash
# run blupf90 programs to estimate gebvs
# running in shell script b/c output files
# are written to working directory with no
# option to change
# and don't want to move R working directory around
# 1 is working directory

# load programs (ceres)
module load blupf90

# move to temp directory
cd $1

# renumber
renumf90 renum.txt > renumbering.out
# estimate gebvs and variance components
airemlf90 renf90.par > estimateGEBVs.out
