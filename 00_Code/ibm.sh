#!/bin/bash -l
#-----------------------------
# Performs the IBM Runs of a LAMMPS Trajectory.
#-----------------------------

#SBATCH --job-name="ibm" ## name that will show up in the queue
#SBATCH --output="01_Logs/%x-%j.out" ## filename of the output; default is slurm-[jobID].out
#SBATCH --error="01_Logs/%x-%j.out" ## file for stderr output
#SBATCH --partition=long
#SBATCH --ntasks=1 ## number of tasks (analyses) to run
#SBATCH --nodes=1
#SBATCH --mail-user s_kell14@uni-muenster.de
#SBATCH --exclude=kaa-1

#SBATCH --mail-type FAIL ## slurm will email you when your job fails
#SBATCH --kill-on-invalid-dep=yes

module purge
module use /home/s_kell14/easybuild/software/modules/all/
module load MYLAMMPS
module load Python/3.11.3-GCCcore-12.3.0

source /mnt/cephflash/s_kell14/venv/ibm/bin/activate

index=$1
nr_steps=$2
temp=$3

l_start=$4
l_end=$5

for f in ./04_Restart/${index}_*.restart; do
    step=${f##*/}
    step=${step#${index}_}
    step=${step%.restart}
    echo $step
done | sort -n > ./${index}_restart_list.txt


for level in $(seq "${l_start}" "${l_end}")
do
    input=./02_PEL/${index}_IS.dat
    python MB_Analysis.py --mode IBM --input "${input}" \
    | while read -r ti tj; do

        tk=$(awk -v ti="$ti" '$1 <= ti' ${index}_restart_list.txt | tail -1)

        mpirun lmp -in con.lmp -var idx ${index} -var nr_steps $(((tj-ti)/2+(ti-tk))) -var num ${tk} -var tempval ${temp} < /dev/null 
    done
done
