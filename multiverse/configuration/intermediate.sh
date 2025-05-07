#!/bin/bash
module load apptainer

container=~/multiverse.sif
batch_sh=~/multiverse/configuration/batch.sh
reproducibility=$2/processed/reproducibility


if [ ! -f $container ]; then
    singularity build multiverse.sif docker://gseasons/multiverse:cluster \
    #singularity pull library://gseasons/analysis/multiverse \
    || echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif"
    exit
fi

singularity exec -B $2:/scratch -B ~/multiverse:/code/multiverse -B $1:/data $container /bin/bash -c "source activate multiverse; export USER=$USER ; python /code/multiverse/run_multiverse.py ${3}"

for filename in $reproducibility/*_workflow_*.pkl; do
    [[ $filename =~ ^.*/(.*)_workflow_(.*).pkl ]] && task=${BASH_REMATCH[1]} && batch=${BASH_REMATCH[2]}
    sbatch --nodes=$4 --ntasks=$5 --account=$6 --time=$7 --mem=$8 $batch_sh $1 $2 $batch $task
    sleep 10
done
