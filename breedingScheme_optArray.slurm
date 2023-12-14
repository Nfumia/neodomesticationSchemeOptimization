#!/bin/bash
#SBATCH --job-name=breedingScheme_opt # Name for your job
#SBATCH --account=macnut_genome
#SBATCH --partition=macnut_genome

#SBATCH --time=7-00:00:00 # Runtime in minutes. Default is 10 minutes. The Maximum runtime currently is 72 hours, 4320 minutes -requests over that time will not run

#SBATCH --nodes=1 # Number of nodes to spread cores across - default is 1 - if you are not using MPI this should likelybe 1
#SBATCH --tasks-per-node=1 # Number of tasks when using MPI. Default is 1
#SBATCH --cpus-per-task=1 # Number of cores requested, Default is 1 (total cores requested = tasks x cores)
#SBATCH --mem=33G # max amount of memory per node you require in MB (default) G (gigabyte)
#SBATCH --array=90-1728

#SBATCH -e /home/nfumia/%x-%A_%a.err # Standard err goes to this file
#SBATCH --mail-user nfumia@hawaii.edu # this is the email you wish to be notified at
#SBATCH --mail-type ALL # this specifies what events you should get an email about ALL will alert you of jobbeginning,completion, failure etc

module purge

parameters=$(sed -n ${SLURM_ARRAY_TASK_ID}p parameters.txt)

module load lang/R/4.2.1-foss-2022a
Rscript breedingSchemeOpt_neoDomestication_1POP_Array.R $parameters