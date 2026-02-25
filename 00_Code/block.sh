#!/bin/bash -l
#-----------------------------
# Runs the later LAMMPS Simulations for the MB Analysis.
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
module load Python/3.11.3-GCCcore-12.3.0

source /mnt/cephflash/s_kell14/venv/ibm/bin/activate

index=$1
nr_steps=$2
temp=$3
start=$4

if python MB_Analysis.py --mode IBM --input "./02_PEL/${index}_IS.dat" | grep -q . ; then
    echo "IBM not finished yet."
    exit 1
fi

latest="./04_Restart/${index}_${start}.restart"
archive="./05_Archive/"
mkdir -p "$archive"

shopt -s nullglob
restart_files=(./04_Restart/${index}_*.restart)
to_archive=()

for f in ${restart_files[@]}; do
    if [[ "$f" != "$latest" ]]; then
        to_archive+=("$f")
    fi
done

if (( ${#to_archive[@]} > 0 )); then
    archive_name="${archive}${index}_block_${start}.tar.gz"
    tar -czf "$archive_name" "${to_archive[@]}"
    rm "${to_archive[@]}"
fi
python MB_Analysis.py --mode Clean --input "./02_PEL/${index}_IS.dat"

y=$((nr_steps/2500))

mpirun lmp -in block.lmp -var idx ${index} -var num ${start} -var nr_steps ${nr_steps} -var tempval ${temp} -var do_restart 1

# When nothing happens anyway, you can also just compare the last and first restart before the IBM:
for (( i=1; i<y; i++)); do
    mpirun lmp -in con.lmp -var idx ${index} -var nr_steps 0 -var num $((i*2500+start)) -var tempval ${temp} -var do_dump 1;
done
