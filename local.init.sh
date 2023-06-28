# initialize env in taiwania for workflow

# set -o pipefail
# set -x
# umask 002

# modules
module purge
module load biology/git/2.34.1
module load libs/singularity/3.10.2

# conda deactivate
# module load biology/snakemake/snakemake
# module load biology/miniconda/miniconda3
conda activate pb-human-wgs
