# AnalysisMultiverse
Software for automated multiverse analysis for fMRI

# Run Instructions
- Launch Terminal
  1. For local run 'python /PATH/TO/MULTIVERSE.PY/multiverse.py -r -d /PATH/TO/BIDS_DATA -o /PATH/TO/OUTPUT_FOLDER'
  2. For cluster run 'python /PATH/TO/MULTIVERSE.PY/multiverse.py -c'
    a. Configure desired settings and parameter space for the cluster (select SLURM for running mode)
    b. After configuration copy files and directories (excluding updated, hide) to compute cluster
    c. Run step 1
  3. To rerun an analysis using the same pipelines as a prior run: 'python /PATH/TO/MULTIVERSE.PY/multiverse.py -r -rr -d /PATH/TO/BIDS_DATA -o /PATH/TO/OUTPUT_FOLDER'
    a. This requires the output directory to contain the reproducibility directory (including generation files) of the analysis to reproduce
    b. Similarly, the files multiverse/configuration/general_configuration.pkl and multiverse/configuration/multiverse_configuration.pkl should be the same as the analysis to be reproduced
  

# Current Run Instructions (for Erin)
- Docker
  1. Launch docker
  2. In terminal run the following: 
    *docker build --tag multiverse --file {PATH_TO_DOCKER_FILE_IN_CODE_DIRECTORY}/multiverse.Dockerfile .*
  3. If running on a Mac, it will be necessary to go to SETTINGS -> RESOURCES and change Memory to at least 8GB (I set mine to 12 to be safe) (NOTE TO GRAHAM: WHAT SETTINGS?)
  4. To launch a container once built:  
    *docker run --rm -it -v {PATH_TO_CODE_DICRECTORY}:/code -v {PATH_TO_BIDS_DATA}:/data -v {PATH_TO_STORE_OUTPUTS}:/scratch -w /scratch multiverse*
    
- Multiverse
  1. Place *expand.pkl*, *master.pkl*, and *population.pkl* in the {PATH_TO_STORE_OUTPUTS} directory
  2. After launching the container, run the following command in the container's terminal: 
    *python /code/run_multiverse.py*
  
  
