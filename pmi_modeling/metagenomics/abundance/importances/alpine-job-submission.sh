#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=MAG_import
#SBATCH --output=alpinejob.MAG-importances.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zachary.burcham@colostate.edu

## make sure you 'module load slurm/alpine' before running job
# purge all existing modules
module purge

# load any modules needed to run your program  
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
source activate ml_pmi

# The directory where you want the job to run
cd /projects/zburcham@colostate.edu/machine_learning/mags/importances

# Job
python MAG_species_importances.py
