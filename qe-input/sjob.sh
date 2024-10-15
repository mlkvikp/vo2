#!/bin/bash
#SBATCH -n 36 -C ib
#SBATCH --time=1:00:00
#SBATCH --job-name=test
#SBATCH --output=std_euler_%j.out
#SBATCH --error=std_euler_%j.err

seedname="m1"
QEDIR="/cluster/project/spaldin/quantum_espresso/7.3/"

mpirun ${QEDIR}pw.x -procs 36 -i $seedname.scf.in > $seedname.scf.out
mpirun ${QEDIR}pw.x -procs 36 -i $seedname.nscf.in > $seedname.nscf.out
mpirun ${QEDIR}projwfc.x -procs 36 -i $seedname.pdos.in > $seedname.pdos.out
mpirun ${QEDIR}pw.x -procs 36 -i $seedname.bands.in > $seedname.bands.out
mpirun ${QEDIR}bands.x -procs 36 -i $seedname.bands.bands.in > $seedname.bands.bands.out
