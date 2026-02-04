#!/bin/bash

echo "WE ARE IN THE BATCH.SH #Mah"
# module load ipykernel/2025a  #Mah - to run the code without the hostenv_multiverse i need to load this module first
module load apptainer #Mah
# module load apptainer/1.2.4 #Mah

container=~/multiverse.sif
custom_base=/opt/miniconda-latest/envs/multiverse/lib/python3.8/site-packages/nipype/pipeline/plugins/base.py
#if templateflow breaks
templates=/home/$USER/.cache/templateflow
ipyth=/scratch_dir/.ipython
# ipyth=/home/scratch_dir/.ipython #/home/$USER/.ipython : #this didnt work #Mah
# export IPYTHONDIR=/home/scratch_dir/.ipython #Mah

###Gra - I got a syntax error here
# if [ ! -f $container ]; then
#     singularity build multiverse.sif docker://gseasons/multiverse:cluster \
#     #singularity pull library://gseasons/analysis/multiverse \
#     || echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif"
#     exit
# fi


###Mah - corrected syntax
if [ ! -f $container ]; then
    #singularity pull library://gseasons/analysis/multiverse \
    singularity build multiverse.sif docker://gseasons/multiverse:cluster || \
    echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif"
    exit
fi

profile=job_${SLURM_JOB_ID}_$(hostname)


echo "Creating profile ${profile}"
ipython profile create ${profile}

echo "Launching controller"
ipcontroller --ip="*" --profile=${profile} --log-to-file &
sleep 45
# sleep 10 #Gra

echo "Launching engines"
#if templateflow breaks
srun singularity run -B $2:/scratch_dir -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/multiverse/templateflow:$templates -B ~/multiverse:/code/multiverse -B $1:/data $container ipengine --profile=${profile} --location=$(hostname) --log-to-file &
#normal
#srun singularity run -B $2:/scratch_dir -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/.ipython:/scratch_dir/.ipython -B ~/multiverse:/code/multiverse -B $1:/data $container ipengine --profile=${profile} --location=$(hostname) --log-to-file &
sleep 45
echo "Launching Job"

### Mah - jun5,2025 ----------------------------------------------------------------------
# echo " MAH| Checking if ipcontroller is running:"
# ps aux | grep ipcontroller

### --------------------------------------------------------------------------------------


#normal
#singularity exec -H $2:/scratch_dir -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/.ipython:/scratch_dir/.ipython -B ~/multiverse:/code/multiverse -B $1:/data \
#if templateflow breaks


singularity exec -B $2:/scratch_dir -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/multiverse/templateflow:/scratch_dir/.cache/templateflow -B ~/multiverse:/code/multiverse -B $1:/data \
$container /bin/bash -c \
"source activate multiverse ; export USER=$USER ; python /code/multiverse/batch_multiverse.py ${3} ${4} ${profile}"

JOB_CHECK=$(sq)

#normal
#singularity exec -H $2:/scratch_dir -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/.ipython:/scratch_dir/.ipython -B ~/multiverse:/code/multiverse -B $1:/data \
#if templateflow breaks


singularity exec -B $2:/scratch_dir -e -B ~/multiverse/plugins_base.py:$custom_base -B ~/multiverse/templateflow:/scratch_dir/.cache/templateflow -B ~/multiverse:/code/multiverse -B $1:/data \
$container /bin/bash -c \
"source activate multiverse ; export USER=$USER ; python /code/multiverse/batch_multiverse_processing.py ${4} '${JOB_CHECK}'"
