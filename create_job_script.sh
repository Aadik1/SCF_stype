#!/bin/bash -l

# --- Slurm Headers ---
#SBATCH --job-name=S2+G5L_FullNode
#SBATCH --account=k19020583
#SBATCH --time=12:00:00
#SBATCH --nodes=1                
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=128      
#SBATCH --mem=0                  
#SBATCH --output=%x.%j.out      

# --- Modules ---
module purge
module load gcc/13.2.0
module load openblas/0.3.24-gcc-13.2.0-multithreading-openmp

# --- Environment Variables ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- Execution ---
/scratch/users/k19020583/

