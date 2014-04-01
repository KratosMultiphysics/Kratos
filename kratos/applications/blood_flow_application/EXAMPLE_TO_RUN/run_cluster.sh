#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test%j.out
#SBATCH --error=test_%j.err
#SBATCH --ntasks=12
#SBATCH --exclusive
#SBATCH --partition=XeonE5645
#SBATCH --nodelist=pez014
#SBATCH --reservation=esoudah_66
 
##Optional - Required memory in MB per node, or per core. Defaults are 3GB per core and 16GB per node.
##SBATCH --mem=2000
##SBATCH --mem-per-cpu=1000
 
##Optional - Estimated execution time
##Acceptable time formats include  "minutes",   "minutes:seconds",
##"hours:minutes:seconds",   "days-hours",   "days-hours:minutes" ,"days-hours:minutes:seconds".
##SBATCH --time=05:00:00
 
########### Further details -> man sbatch ##########

export OMP_NUM_THREADS=1
export PYTHONPATH=$PYTHONPATH:/home/esoudah/kratos/ #Aqui va el path on hagis compilar cratos

##module load openmpi/1.6.2

##time python script.py
time python ToRun_new_27032014.py
