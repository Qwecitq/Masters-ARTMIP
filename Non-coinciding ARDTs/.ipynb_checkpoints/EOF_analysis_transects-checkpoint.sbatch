#!/bin/bash
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ktquagra@iu.edu
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 03:00:00
#SBATCH -A m1517


bash EOF_analysis_transects.sh