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

index=$1
nr_steps=$2
temp=$3

l_start=$4
l_end=$5

for level in $(seq "${l_start}" "${l_end}")
do
	t0=$(date +%s)
    input=./02_PEL/${index}_IS.dat

    uv run MB_Analysis.py --mode IBM --input "${input}" \
    | while read -r ti tj; do
	
	t1=$(date +%s)

	tk=$(for f in ./04_Restart/${index}_*.restart; do
	        step=${f##*/}
        	step=${step#${index}_}
        	step=${step%.restart}
        echo "$step"
	done | awk -v ti="$ti" '$1 <= ti' | sort -n | tail -1)	
	
	t2=$(date +%s)

	mpirun lmp -in con.lmp -var idx ${index} -var nr_steps $(((tj-ti)/2+(ti-tk))) -var num ${tk} -var tempval ${temp} < /dev/null 
    
	t3=$(date +%s)

	echo "TIMING ti=$ti tj=$tj | shell=$((t2-t1))s | lammps=$((t3-t2))s" >> ./timing.log
	
	done

	t4=$(date +%s)
	echo "TOTAL MB_Analysis + loop: $((t4-t0))" >> ./timing.log

done
