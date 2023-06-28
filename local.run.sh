source /work/hsinhan97/ncku_human/template/workflow/local.init.sh

# step 1
bash workflow/process_smrtcells.slurm.sh ; echo "END" ; sleep 7d

# step 2
bash workflow/process_sample.slurm.sh sample_A1 ; echo "END" ; sleep 7d
bash workflow/process_sample.slurm.sh sample_B1 ; echo "END" ; sleep 7d

# step 3
bash workflow/process_cohort.slurm.sh cohort_A1 ; echo "END" ; sleep 7d
bash workflow/process_cohort.slurm.sh cohort_B1 ; echo "END" ; sleep 7d
