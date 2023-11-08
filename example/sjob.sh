#!/bin/bash

#SBATCH -n 36
#SBATCH --time=01:00:00
#SBATCH --job-name=ac200
#SBATCH --output=std_euler_%j.out
#SBATCH --error=std_euler_%j.err

mpirun pw.x -procs 36 -i m1.scf.in > m1.scf.out
mpirun pw.x -procs 36 -i m1.bands.in > m1.bands.out 
mpirun bands.x -procs 36 -i m1.bands.bands.in > m1.bands.bands.out




