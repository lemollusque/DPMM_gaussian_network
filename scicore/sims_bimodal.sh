#!/bin/bash

#SBATCH --job-name=sims                   #This is the name of your job
#SBATCH --cpus-per-task=16                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=1G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=out_bimodal.txt     #These are the STDOUT and STDERR files
#SBATCH --error=err_bimodal.txt
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=bastian.widmer@unibas.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
module load R/4.2.1-foss-2022a 

#export your required environment variables below
#################################################


#add your command lines below
#############################
Rscript sims_bimodal.R
