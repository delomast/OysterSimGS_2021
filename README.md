# OysterSimGS_2021
Simulations to compare genotyping strategies for genomic selection in 
oyster breeding programs

This repository holds the code used to perform the simulations described in _to be updated_.

generate_random_seeds.R is an R script that generates the random seeds used by the simulations. 
This was run once prior to running the simulations. The simulations were run on an HPC 
using the Slurm Workload Manager as array jobs (one task for each iteration). The following is an explanation of which files 
correspond to each scenario. 

Base scenario

- scrm_SNPsim_arrary.sh is a shell script that launches the simulation
- multGen_scrm_ai2_HPC.R is an R script that executes the simulation and saves the results

Random SNP selection

- scrm_SNPsim_arrary_rand.sh is a shell script that launches the simulation
- multGen_scrm_ai2_HPC_rand.R is an R script that executes the simulation and saves the results

Small family size

- scrm_SNPsim_arrary_10off.sh is a shell script that launches the simulation
- multGen_scrm_ai2_HPC_numOff.R is an R script that executes the simulation and saves the results

Current generation only

- curGen_SNPsim_arrary.sh is a shell script that launches the simulation
- onlyCurGen_scrm.R is an R script that executes the simulation and saves the results

Empirical founder populations

- snp_empir_HPC.R is an R script that executes the simulation and saves the results
- empirMBP_SNPsim_arrary.sh is a shell script that launches the MBP Pacific oysters simulation
- empirNgulf_SNPsim_arrary.sh is a shell script that launches the Northern Gulf eastern oysters simulation
- empirFRC_SNPsim_arrary.sh is a shell script that launches the French Pacific oysters simulation
- empirENG_SNPsim_arrary.sh is a shell script that launches the British Isles Pacific oysters simulation
- empirCAN_SNPsim_arrary.sh is a shell script that launches the Canadian eastern oysters simulation


The other files in the repository are:

- utils.R : functions used by the simulation scripts
- replaceArrayNames.py : replaced chr and locus names in array data with data from mapping the loci to the reference genome
- runAlphaImpute2.sh : impute with alphaImpute2
- run_blupf90.sh : calculate GEBVs with airemlf90
- sim_oyster_genos_scrm.sh : simulate founder haplotypes with scrm
- results.R : summarizing the outputs, making figures and tables
- archive/* : old files not used in the final simulations



This software is a "United States Government Work" under the terms of the United States Copyright Act. 
It was written as part of the authors' official duties as United States Government employees and 
thus cannot be copyrighted. This software is freely available to the public for use.
