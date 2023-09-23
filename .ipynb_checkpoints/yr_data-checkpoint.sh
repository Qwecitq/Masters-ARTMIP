#!/bin/bash
#SBATCH -C knl
#SBATCH -N 73
#SBATCH -q regular
#SBATCH -t 00:30:00
#SBATCH -A m1517
#SBATCH -J 2_daily_average

# load the gcc environment
module swap PrgEnv-intel PrgEnv-gnu

# bring a TECA install into your environment
module use /global/cscratch1/sd/loring/teca_testing/installs/develop/modulefiles
module load teca

# print the commands as they execute, and error out if any one command fails
set -e
set -x

#data directory
data_dir=/global/project/projectdirs/m1517/cascade/external_datasets/ARTMIP/tier1/ftp_mirror/cascade_bard_v1

# make a directory for the output files
out_dir=/global/u1/k/kquagra/kquagra/yr_files/cascade_bard_v1
mkdir -p ${out_dir}

# compute the daily average. change -N and -n to match the run size.
# the run size is determened by the number of output time steps. here the
# input is 3 hourly, the output is daily.
time srun -N 73 -n 146 \
    teca_temporal_reduction \
        --n_threads 2 --verbose 1 --input_regex ${data_dir}/'.*\.nc$' \
        --interval monthly --operator average --point_arrays ar_binary_tag --ignore_fill_value \
        --output_file ${out_dir}/MERRA_cascade_bard_monthly%t%.nc \
        --file_layout yearly