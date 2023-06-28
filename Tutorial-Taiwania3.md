# How to run this on Taiwania 3

# Preface
## Taiwania 3
- Managed by slurm system.
- Heavy computing job is not allowed on login node.
- Get kicked out if session is idle for a couple minutes. Use the magic `sleep 7d`.
- Built-in git is so outdated. If you got mysterious errors during installing something, try `module load biology/git/2.34.1`.
- 100 concurrent jobs per user at maximum.
- Do not support `ftp` or `lftp`, try `rclone`.
- Use `/work/username` for computing instead of `/home/username`.

## Snakemake
- Was designed for user to submit a master job, and the master job will submit resting computing jobs.
- Some light computing jobs may run aside master job, so called `localrules`.
- Because Taiwania do not allow submitting job from computing node. Here we are going to run master job on login node, and let it submit all computing jobs to computing node.
- All `localrules` are disabled to fit Taiwania policy.

# Setup
## Workspace
- Prepare working directory structure according to https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/Tutorial.md#2-prepare-workspace.
- `reference/` and `resources/` are provided by PacBio.
```sh
DIR=/path/to/somewhere/

# clone workflow, checkout `twnia3` branch
cd $DIR/template
git clone -b twnia3 --single-branch https://github.com/cgjosephlee/pb-human-wgs-workflow-snakemake

# create working directories, edit script before run
cd $DIR/work
mkdir run01
cd run01
bash $DIR/template/pb-human-wgs-workflow-snakemake/local.setup.sh

# create sample folder and link raw sequencing files
mkdir -p smrtcells/ready/sample_A1
ln -s /path/to/xxx_reads.bam smrtcells/ready/sample_A1
```

## Configurations
- https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/Tutorial.md#3-analysis-configuration
- Edit `config.yml` if necessary.
- Create `cohort.yml`. If singleton, every sample is a distinct cohort.

## Environment
- Use conda env. Both user installed or use conda module in Taiwania are valid.
- Git is require for installing dependencies. Build-in git version is obsoleted, use `module load biology/git/2.34.1`.
- Conda env for computing jobs are created automatically by snakemake.
- Require singularity.
```sh
# create env for master job
module load biology/git/2.34.1
module load biology/miniconda/miniconda3  # use conda module
mamba create -c conda-forge -c bioconda -n pb-human-wgs snakemake=6.15.3 tabulate=0.8.10 pysam=0.16.0.1 python=3

# initialize local env
bash workflow/local.init.sh
```

# Run
- What master job do is to monitor the submitted jobs and output files only, no computing is run on login node so it won't be killed.
- Use `tmux`, `screen` or `nohup` to run snakemake master job. I used `tmux`.
- Commands are listed in `local.run.sh`. They should be run step by step.
- Ended by `sleep 7d` to halt the session so we can check logs when job is finished.
```sh
# step 1
bash workflow/process_smrtcells.slurm.sh ; echo "END" ; sleep 7d

# step 2
# run on each sample, can be run in parallel
bash workflow/process_sample.slurm.sh sample_A1 ; echo "END" ; sleep 7d
bash workflow/process_sample.slurm.sh sample_B1 ; echo "END" ; sleep 7d

# step 3
# run on each cohort, can be run in parallel
bash workflow/process_cohort.slurm.sh cohort_A1 ; echo "END" ; sleep 7d
bash workflow/process_cohort.slurm.sh cohort_B1 ; echo "END" ; sleep 7d
```

# Outputs
- Check https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake/blob/main/Tutorial.md#outputs.
