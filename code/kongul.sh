!#/bin/bash

# The PBS lines are for the que system.
#PBS -N test_kongull
# "default" refers to the kind of CPU. alternative is "intel". The ppn refers to _MPI_ processes per node.
#PBS -lnodes=4:ppn=12:default
# Max time. Kill afterwards
#PBS -lwalltime=00:00:30
# Memory per MPI process. 24 GB per node available.
#PBS -lpmem=2000MB
# Acount (who pays)
#PBS -A freecycle
# que identity.
#PBS -q optimist
# what to do with output. (stdout and stderr in same file)
#PBS -j oe

# This is a bash script (run once on one node)
cd ${PBS_O_WORKDIR}

# To get access to libraries. Demands the compilators to be loaded.
module load intelcomp
module load openmpi/1.4.3-intel
# How memory is used during multithreading.
KMP_AFFINITY="granularity=fine , compact"

for n in `seq 3 14`
do
	ELAPSED=`OMP_NUM_THREADS=1 mpirun -npernode 1 ./oppg1 1 $n | grep "^[0-9]"`
 	echo "$ELAPSED" >> results.mat
done

