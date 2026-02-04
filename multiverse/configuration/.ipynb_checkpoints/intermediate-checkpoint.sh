#!/bin/bash

module load apptainer #Mah
# module load apptainer/1.2.4 #Mah

container=~/multiverse.sif
batch_sh=~/multiverse/configuration/batch.sh
reproducibility=$2/processed/reproducibility
echo "we are in intermediate.sh |Mah "
###Gra # I get a syntax error here
# if [ ! -f $container ]; then
#     singularity build multiverse.sif docker://gseasons/multiverse:cluster \
#     #singularity pull library://gseasons/analysis/multiverse \
#     || echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif"
#     exit
# fi


###MAh # corrected version
if [ ! -f $container ]; then
    # singularity pull library://gseasons/analysis/multiverse
    singularity build multiverse.sif docker://gseasons/multiverse:cluster || \
    echo "Cannot access container, please upload the image into your home directory: https://cloud.sylabs.io/library/gseasons/multiverse/multiverse.sif"
    exit
fi

export REQUESTS_CA_BUNDLE=$PWD/cacert.pem #|Mah

singularity exec -B $2:/scratch_dir -B ~/multiverse:/code/multiverse -B $1:/data $container /bin/bash -c "source activate multiverse; export USER=$USER ; python /code/multiverse/run_multiverse.py ${3}"

#### original run! ######
# for filename in $reproducibility/*_workflow_*.pkl; do
#     [[ $filename =~ ^.*/(.*)_workflow_(.*).pkl ]] && task=${BASH_REMATCH[1]} && batch=${BASH_REMATCH[2]}
#     sbatch --nodes=$4 --ntasks=$5 --account=$6 --time=$7  --mem=$8 --mail-user=mahshid.soleymani@ucalgary.ca --mail-type=ALL $batch_sh $1 $2 $batch $task
#     # sleep 10 #Gra
#     sleep 45 #Mah
# done
##########################

#### run1 : run only one pipeline; rest_workflow_8 : batch#17 #######
for filename in $reproducibility/rest_workflow_15.pkl; do
    [[ $filename =~ ^.*/(.*)_workflow_(.*).pkl ]] && task=${BASH_REMATCH[1]} && batch=${BASH_REMATCH[2]}
    sbatch --nodes=$4 --ntasks=$5 --account=$6 --time=20:10:00  --mem=$8 --mail-user=mahshid.soleymani@ucalgary.ca --mail-type=ALL $batch_sh $1 $2 $batch $task
    # sleep 10 #Gra
    sleep 45 #Mah
done
#########################################

echo $1 "1"
echo $2 "2"
echo $3 "3"
echo $4 "4"
echo $5 "5"
echo $6 "6"
echo $7 "7"
echo $8 "8"
echo $batch "batch"
echo $task "task"
echo ${BASH_REMATCH[1]}  "{BASH_REMATCH[1]} or taks"
echo ${BASH_REMATCH[2]}  "{BASH_REMATCH[2]} or batch"

### --time=04:00:00 (for testing)
### change it back to --time=$7 , --mem=$8 , --ntasks=$5

