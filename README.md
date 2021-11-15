# AnalysisMultiverse
Software for automated multiverse analysis for fMRI

# Current Run Instructions (for Erin)
- Docker
  1. Launch docker
  2. In terminal run the following: 
    *docker build --tag multiverse --file {PATH_TO_DOCKER_FILE_IN_CODE_DIRECTORY}/multiverse.Dockerfile .*
  3. To launch a container once built:  
    *docker run --rm -it -v {PATH_TO_CODE_DICRECTORY}:/code -v {PATH_TO_BIDS_DATA}:/data -v {PATH_TO_STORE_OUTPUTS}:/scratch -w /scratch multiverse*
    
- Multiverse
  1. After launching the container, run the following command in the container's terminal: 
    *python /code/run_multiverse.py*
  
  
