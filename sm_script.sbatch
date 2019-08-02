#!/bin/sh

# sometimes if something is running for a while the powers that be kill it, so start a interactive job...
# sdev -p normal -t 34:59:59 # should check in and kill before time runs out

#SBATCH --job-name=sm_master
#SBATCH --output=SM_%j_pilot.out
#SBATCH --error=SM_%j_pilot.err
#SBATCH --time=30:00:00
#SBATCH -p pritch,normal,owners
#SBATCH --mem=15G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
###SBATCH --mail-user=antmarge@stanford.edu
	#=aow run normal batch commands

#snakemake --unlock
#20170630/3.bwa/bwa.RMPR-11-U_1-50580648.bam
# took out owners for now so doesn't get reset

snakemake --rerun-incomplete --timestamp --keep-going -j 500 --cluster-config sm_slurm_config.json --cluster "sbatch -p {cluster.queue} --ntasks-per-node={cluster.tasks} --job-name={cluster.name} -o {cluster.error} -e {cluster.error} -t {cluster.time} --nodes={cluster.nodes} --mem={cluster.memory}"
