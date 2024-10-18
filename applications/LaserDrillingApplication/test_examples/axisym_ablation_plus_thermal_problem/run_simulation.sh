#!/bin/bash
#SBATCH --job-name=axisym_ablation_plus_thermal_problem
#SBATCH --output=axisym_ablation_plus_thermal_problem.out
#SBATCH --error=axisym_ablation_plus_thermal_problem.err
#SBATCH --partition=Piksel
#SBATCH --ntasks-per-node=20

##Optional - Required memory in MB per node, or per core. Defaults are 1GB per core.
##SBATCH --mem=3096
#SBATCH --mem-per-cpu=3096

##Optional - Estimated execution time
##Acceptable time formats include  "minutes",   "minutes:seconds",
##"hours:minutes:seconds",   "days-hours",   "days-hours:minutes" ,"days-hours:minutes:seconds".
#SBATCH --time=10-0

########### Further details -> man sbatch ##########

export PYTHONPATH=/home/slatorre/lasers/Kratos/bin/Release/:$PYTHONPATH
export LD_LIBRARY_PATH=/home/slatorre/lasers/Kratos/bin/Release/libs:$LD_LIBRARY_PATH
export PATH=/opt/intel/oneapi/mkl/2022.2.1/lib/intel64:$PATH

export OMP_NUM_THREADS=20
/opt/intel/oneapi/intelpython/python3.9/bin/python3 MainKratos.py

