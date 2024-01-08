#PBS -N edison_single_cell_test
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

rm -rf *.plt *.png result/*.plt

./cell -input_deck input_single.txt > logfile.txt

