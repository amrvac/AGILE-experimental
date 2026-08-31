#!/bin/bash
#SBATCH --partition=gpu_h100
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH -t 00:10:00
##SBATCH --mem=0
#SBATCH --job-name=check
#SBATCH -o slurms/single_out-%j # STDOUT 
##SBATCH -e slurms/single_err-%j # STDERR

module purge
module load 2025
module load OpenMPI/5.0.7-NVHPC-25.3-CUDA-12.8.0
source $AGILE_DIR/.venv/bin/activate

make clean-all
make -j arch=nvidia OPENACC=1 USE_MPIWRAPPERS=1

srun ./agile -i agile.par

###SBATCH -N 1
###SBATCH --gpus-per-node 1
###SBATCH --cpus-per-task 1
###SBATCH --mem=250G
###SBATCH --exclusive
####SBATCH --gpus-per-node 2
####SBATCH --ntasks-per-node 2
####SBATCH --gpu-bind=map_gpu:0,1,2,3
###SBATCH -t 00:12:00
###SBATCH -o out
