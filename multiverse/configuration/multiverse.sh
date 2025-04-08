#!/bin/bash

module load apptainer

container=~/multiverse.sif
custom_base=/opt/miniconda-latest/envs/multiverse/lib/python3.8/site-packages/nipype/pipeline/plugins/base.py
#if templateflow breaks
templates=/home/$USER/.cache/templateflow

if [ ! -f $container ]; then
    singularity build multiverse.sif docker://gseasons/multiverse:cluster \
    #singularity pull library://gseasons/analysis/multiverse \
    || echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/analysis/multiverse.sif"
    exit
fi

profile=job_${SLURM_JOB_ID}_$(hostname)

echo "Creating profile ${profile}"
ipython profile create ${profile}

echo "Launching controller"
ipcontroller --ip="*" --profile=${profile} --log-to-file &
sleep 10

echo "Launching engines"
#if templateflow breaks
srun singularity run -B $2:/scratch -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/multiverse/templateflow:$templates -B ~/.ipython:/scratch/.ipython -B ~/multiverse:/code/multiverse -B $1:/data $container ipengine --profile=${profile} --location=$(hostname) --log-to-file &
#normal
#srun singularity run -B $2:/scratch -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/.ipython:/scratch/.ipython -B ~/multiverse:/code/multiverse -B $1:/data $container ipengine --profile=${profile} --location=$(hostname) --log-to-file &
sleep 45
echo "Launching Job"
#normal
#singularity exec -H $2:/scratch -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/.ipython:/scratch/.ipython -B ~/multiverse:/code/multiverse -B $1:/data \
#if templateflow breaks
singularity exec -H $2:/scratch -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/multiverse/templateflow:/scratch/.cache/templateflow -B ~/.ipython:/scratch/.ipython -B ~/multiverse:/code/multiverse -B $1:/data \
$container /bin/bash -c \
"source activate multiverse ; export USER=$USER ; python /code/multiverse/run_multiverse.py ${3} ${profile}"
