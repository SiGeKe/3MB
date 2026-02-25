#!/bin/bash -l
#-----------------------------
# Performs the Analysis on an IS Trajectory.
#-----------------------------

#SBATCH --job-name="analysis" ## name that will show up in the queue
#SBATCH --output="01_Logs/%x-%j.out" ## filename of the output; default is slurm-[jobID].out
#SBATCH --error="01_Logs/%x-%j.out" ## file for stderr output
#SBATCH --partition=short
#SBATCH --ntasks=1 ## number of tasks (analyses) to run
#SBATCH --nodes=1
#SBATCH --mail-user s_kell14@uni-muenster.de
#SBATCH --exclude=kaa-1

#SBATCH --mail-type FAIL ## slurm will email you when your job fails
#SBATCH --kill-on-invalid-dep=yes

module purge
module use /home/s_kell14/easybuild/software/modules/all/
module load Python/3.11.3-GCCcore-12.3.0

source /mnt/cephflash/s_kell14/venv/ibm/bin/activate

archive="./06_Results/"
mkdir -p "$archive"

python MB_Analysis.py --input "./02_PEL/*_IS.dat" --mode Lifetimes --output "./06_Results/Lifetimes.dat"

for file in ./02_PEL/*_IS.dat; do
    base=$(basename "$file" _IS.dat)
    python MB_Analysis.py --input "$file" --mode MB --output "./06_Results/${base}_MB.dat"
done
