#!/bin/bash

module load singularity

container=~/multiverse_latest.sif
custom_base=/opt/miniconda-latest/envs/multiverse/lib/python3.8/site-packages/nipype/pipeline/plugins/base.py

if [ ! -f $container ]; then
    singularity pull library://gseasons/analysis/multiverse \
    || echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif"
    exit
fi

profile=job_${SLURM_JOB_ID}_$(hostname)

echo "Creating profile ${profile}"
ipython profile create ${profile}

echo "Launching controller"
ipcontroller --ip="*" --profile=${profile} --log-to-file &
sleep 10

echo "Launching engines"
srun singularity run -B $2:/scratch -e -B ~/code/plugins_base.py:$custom_base -B ~/.ipython:/scratch/.ipython -B ~/code:/code/multiverse -B $1:/data $container ipengine --profile=${profile} --location=$(hostname) --log-to-file &
sleep 45
echo "Launching Job"

singularity exec -H $2:/scratch -e -B ~/code/plugins_base.py:$custom_base -B ~/.ipython:/scratch/.ipython -B ~/code:/code/multiverse -B $1:/data \
$container /bin/bash -c \
"source activate multiverse ; python /code/multiverse/run_multiverse.py ${3} ${profile}"
