#!/bin/bash

#SBATCH -n 36
#SBATCH --time=1:00:00
#SBATCH --job-name=m1
#SBATCH --output=std_euler_%j.out
#SBATCH --error=std_euler_%j.err

seedname="m1"

mpirun pw.x -procs 36 -i $seedname.scf.in > $seedname.scf.out
mpirun pw.x -procs 36 -i $seedname.nscf.in > $seedname.nscf.out
mpirun projwfc.x -procs 36 -i $seedname.pdos.in > $seedname.pdos.out
mpirun pw.x -procs 36 -i $seedname.bands.in > $seedname.bands.out
mpirun bands.x -procs 36 -i $seedname.bands.bands.in > $seedname.bands.bands.out
