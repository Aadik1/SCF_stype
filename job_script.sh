#!/bin/bash -l
#$ -N SpinJunction
#$ -P Gold
#$ -A KCL_Kantorovitch
#$ -l h_rt=12:00:00
#$ -l mem=4.5G
#$ -pe smp 40
#$ -cwd

module unload gcc-libs
module load gcc-libs/10.2.0
module load openblas/0.3.13-openmp/gnu-10.2.0

export OMP_NUM_THREADS=40

/home/mmm1550/SCF_stype/run.out
