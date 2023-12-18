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

snakemake --cluster "sbatch -A {cluster.account} -q {cluster.qos} -c {cluster.cpus-per-task} -N {cluster.Nodes} \
  -t {cluster.runtime} --mem {cluster.mem} -J {cluster.jobname} --mail-type={cluster.mail_type} \
  --mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" \
  --cluster-config config/cluster.json --jobs 100 --latency-wait 20 --rerun-incomplete --use-conda