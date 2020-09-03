
#!/bin/bash
####PBS -q global
#PBS -l nodes=m05:ppn=48
#PBS -N Ef202off

cd $PBS_O_WORKDIR
date


OUTFILE=Ef202off.out
export OMP_NUM_THREADS=1
export PYTHONPATH="${PYTHONPATH}$HOME/.local/lib/python3/site-packages"


time mpirun -n 48 python3 1dwire_current_Ef202_noqpc.py > $OUTFILE

date

