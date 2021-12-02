# AnalysisMultiverse
Software for automated multiverse analysis for fMRI

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
  
  
