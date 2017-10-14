#!/bin/bash

#SBATCH --time=360

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=west

#SBATCH --output=nbody.out
#SBATCH --error=nbody.err

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" = 'abu' ]
then
        export MPICH_NEMESIS_NETMOD=ib
fi

srun hostname
time mpiexec ./nbody 1000 0.01 1
