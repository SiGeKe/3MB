#!/bin/bash

name="Test"

index=$1
nr_steps=40000000
temp=0.1333

block=0 # Remember to change!

# 1. Run-Job

# 1.1 For the first block
jid_run=$(sbatch --job-name="${name}-run-${index}-T${temp}" run.sh ${index} ${nr_steps} ${temp} | awk '{print $4}')

# 1.2 For the subsequent blocks
#jid_run=$(sbatch --job-name="${name}-block-${block}-run-${index}-T${temp}" block.sh ${index} ${nr_steps} ${temp} $((nr_steps*block)) | awk '{print $4}')

# 2. Loop of IBM-Jobs

# 2.1 When Starting the IBM
prev_jid=$(sbatch --job-name="${name}-ibm-${index}-T${temp}" --dependency=afterok:${jid_run} ibm.sh ${index} ${nr_steps} ${temp} 1 5 1 | awk '{print $4}')

# 2.2 When Continuuing the IBM
#prev_jid=$(sbatch --job-name="${name}-ibm-${index}-T${temp}" ibm.sh ${index} ${nr_steps} ${temp} 1 5 0 | awk '{print $4}')

for start in 6 11 16; do
    end=$((start+4))
    jid_ibm=$(sbatch --job-name="${name}-ibm-${index}-T${temp}" --dependency=afterok:${prev_jid} ibm.sh ${index} ${nr_steps} ${temp} ${start} ${end} 0 | awk '{print $4}')
    prev_jid=${jid_ibm}
done
