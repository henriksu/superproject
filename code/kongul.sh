!#/bin/bash

#PBS -N test_kongull
#PBS -lwalltime=05:00:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine , compact"

for n in `seq 3 14`
do
	ELAPSED=`OMP_NUM_THREADS=4 mpirun -npernode 12 ./oppg1 $n | grep error | awk -F ' ' '{print $5}'`
	let N=n
 	echo "$N $ELAPSED" >> results.mat
done

