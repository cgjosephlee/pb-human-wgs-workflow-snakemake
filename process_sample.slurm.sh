#!/bin/bash
#SBATCH -A 100humans
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# USAGE: sbatch workflow/process_sample.slurm.sh <sample_id>

SAMPLE=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
SAMPLEDIR="cohorts/${SAMPLE}"
LOCKFILE="$SAMPLEDIR/process_cohort.lock"

if [ -f "$LOCKFILE" ]; then
    echo "lockfile $LOCKFILE already exists. Remove lockfile and try again." && exit 1
else
    mkdir -p "$SAMPLEDIR" || exit 1
    touch "$LOCKFILE" || exit 1
fi
trap 'rm -f ${LOCKFILE}; exit' SIGINT SIGTERM ERR EXIT

# get variables from workflow/variables.env
source workflow/variables.env

# execute snakemake
snakemake \
    --config "sample='$SAMPLE'" \
    --nolock \
    --profile workflow/profiles/slurm \
    --snakefile workflow/process_sample.smk \
    --singularity-args "-B /work"
