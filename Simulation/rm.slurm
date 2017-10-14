#!/bin/bash

#SBATCH --time=360

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=west

#SBATCH --output=tar.out
#SBATCH --error=tar.err

if [ "${SLURM_PARTITION}" = 'abu' ]
then
        export MPICH_NEMESIS_NETMOD=ib
fi

srun hostname
time rm -r <dir>
