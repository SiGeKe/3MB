#!/bin/bash -l
#-----------------------------
# Runs the inital LAMMPS Simulation for the MB Analysis.
#-----------------------------

#SBATCH --job-name="run" ## name that will show up in the queue
#SBATCH --output="01_Logs/%x-%j.out" ## filename of the output; default is slurm-[jobID].out
#SBATCH --error="01_Logs/%x-%j.out" ## file for stderr output
#SBATCH --partition=short
#SBATCH --ntasks=1 ## number of tasks (analyses) to run
#SBATCH --nodes=1
#SBATCH --mail-user s_kell14@uni-muenster.de

#SBATCH --mail-type BEGIN  ## slurm will email you when your job starts
#SBATCH --mail-type END  ## slurm will email you when your job ends
#SBATCH --mail-type FAIL ## slurm will email you when your job fails
#SBATCH --kill-on-invalid-dep=yes

module purge
module use /home/s_kell14/easybuild/software/modules/all/
module load MYLAMMPS

index=$1
nr_steps=$2
temp=$3

y=$((nr_steps/2500))

# When performing simulations starting from MC samples:
mpirun lmp -in run.lmp -var idx ${index} -var nr_steps ${nr_steps} -var tempval ${temp} -var do_restart 1

# When performing simulations starting from LAMMPS files:
# mpirun lmp -in block.lmp -var idx ${index} -var nr_steps ${nr_steps} -var tempval ${temp} -var do_restart 1 -var do_start 1

for (( i=1; i<y; i++ )); do
    mpirun lmp -in con.lmp -var idx ${index} -var nr_steps 0 -var num $((i*2500)) -var tempval ${temp} -var do_dump 1;
done

