#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A cheme                 # The account name for the job.
#SBATCH -J tiffCreate    		 # The job name.
#SBATCH -t 0:10:00               # The time the job will take to run.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb        # The memory the job will use per cpu core.
#SBATCH --array=6-10,26-30

module load matlab/2019b
matlab -nosplash -nodisplay -nodesktop -r "ringplot_prepare_tiff($SLURM_ARRAY_TASK_ID,150)" 

# End of script
