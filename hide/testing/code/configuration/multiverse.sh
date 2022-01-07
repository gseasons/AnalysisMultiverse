#!/bin/bash
#SBATCH --account=def-emazerol
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:00
#SBATCH --mem=8G

cd $2
module load singularity

singularity exec -H $2:/scratch -e -B ~/code:/code/multiverse -B $1:/data \
library://gseasons/multiverse/multiverse.sif:latest /bin/bash -c \
"source activate multiverse ; python /code/multiverse/run_multiverse.py"\
|| singularity exec -H $2:/scratch -e -B ~/code:/code/multiverse -B $1:/data \
~/multiverse.sif_latest.sif /bin/bash -c \
"source activate multiverse ; python /code/multiverse/run_multiverse.py"\
|| singularity exec -H $2:/scratch -e -B ~/code:/code/multiverse -B $1:/data \
multiverse.sif_latest.sif /bin/bash -c \
"source activate multiverse ; python /code/multiverse/run_multiverse.py"\
|| (singularity pull library://gseasons/multiverse/multiverse.sif:latest &&\
singularity exec -H $2:/scratch -e -B ~/code:/code/multiverse -B $1:/data \
multiverse.sif_latest.sif /bin/bash -c \
"source activate multiverse ; python /code/multiverse/run_multiverse.py")\
|| echo "Error: Cannot run or pull singularity container. \
Please upload the container image (https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif)\
 to your working (output) directory or your home directory."