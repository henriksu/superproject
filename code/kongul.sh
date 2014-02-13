#!/bin/bash

# The PBS lines are for the que system.
#PBS -N test_kongull
# "default" refers to the kind of CPU. alternative is "intel". The ppn refers to _MPI_ processes per node.
#PBS -lnodes=4:ppn=12:default
# Max time. Set it as low as possible to get priority. But remember that the job is killed at timeout.
#PBS -lwalltime=00:00:30
# Memory per MPI process. There is 24 GB available on each node.
#PBS -lpmem=2000MB
# Account (who pays)
#PBS -A freecycle
# What que to place it in.
#PBS -q optimist
# what to do with output (stdout and stderr in same file).
#PBS -j oe

# This is a bash script (run once on one node)
cd ${PBS_O_WORKDIR}

# To get access to libraries, the following modules must be loaded.
module load intelcomp
module load openmpi/1.4.3-intel

# How threads are treated. Something about keeping each thread on one processor to avoid copying cache or something like that.
KMP_AFFINITY="granularity=fine , compact"


# From this point on it is an ordinary bash-script, calling the interesting program.
for n in `seq 3 14`
do
	ELAPSED=`OMP_NUM_THREADS=1 mpirun -npernode 1 ./oppg1 1 $n | grep "^[0-9]"`
 	echo "$ELAPSED" >> results.mat
done

