#!/bin/bash
#SBATCH --job-name=TLS-disp
#SBATCH --account=<account>
#SBATCH --qos=<account qos>
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=96:00:00
#SBATCH --output=logs/%j_disp.log
#SBATCH --error=logs/%j_disp.log

# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date


##----------------------------------------------------------
# Modules
module load conda
module load snakemake


##----------------------------------------------------------
# Run

pip show snakemake-executor-plugin-cluster-generic 1>/dev/null
if [ $? == 0 ]; then
   echo "Good to go" #Replace with your actions
else
   pip install snakemake-executor-plugin-cluster-generic
fi

snakemake --profile profiles/slurm/