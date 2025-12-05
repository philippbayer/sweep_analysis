#!/bin/bash --login
#SBATCH --account=pawsey0964
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --partition=work


source /software/projects/pawsey0390/pbayer/conda/etc/profile.d/conda.sh

conda activate ./doWork

srun -c 32 python findSweeps.py vcf_dir 32  > log 2> err
