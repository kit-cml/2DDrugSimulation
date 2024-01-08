#PBS -N single_cell_test
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20000:00:00
#PBS -e stderr.log
#PBS -o stdout.log
# Specific the shell types
#PBS -S /bin/bash
# Specific the queue type
#PBS -q dque


# Script for running the single cell simulation
cd $PBS_O_WORKDIR

ulimit -s unlimited

NPROCS=`wc -l < $PBS_NODEFILE`

echo This job has allocated $NPROCS nodes

cd $PBS_O_WORKDIR

rm -rf *.plt *.png *.txt

#/opt/Lib/openmpi-2.1.3/bin/mpirun -v -machinefile $PBS_NODEFILE -np $NPROCS  
./cell \
    -variant ori \
    -simulation_mode 0 \
    -dt 0.1 \
    -bcl_init 750.0\
    -num_pace1 10 \
    -scaled_channels g_Kr*1.0,g_K1*1.0,g_bca*1.0,g_Na*1.0 \
    -stim_amp 1.0 \
    -stim_dur 1.0 >& logfile.txt

