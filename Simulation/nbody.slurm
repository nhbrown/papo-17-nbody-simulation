#!/bin/bash

#SBATCH --time=360

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=west

#SBATCH --output=nbody.out
#SBATCH --error=nbody.err

if [ "${SLURM_PARTITION}" = 'abu' ]
then
        export MPICH_NEMESIS_NETMOD=ib
fi

srun hostname
time ./nbody 1000 0.01 1
