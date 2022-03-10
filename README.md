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
  4. Note: It may take a long time from the workflow calling .run() to actual execution (connecting nodes/inputs/outputs) - important in multiverse, as the more nodes, subjects, and pipelines the longer this will take (important to ask if some way to request fewer resources while this happening, then step up to needed amount once starts running)
  
# Node Naming Conventions/Multiverse Modification
- To add more parameters to multiverse analysis, edit the default.json file in configuration
  1. Find the specific node in the workflow, cross reference with changeable options for nipype interface
  2. Add new entry in format seen in json file 
    a. Non-numeric parameters need a value_map parameter to indicate it will be aliased in the genetic algorithm
    b. Default category indicates the index of the default value
    c. alias is what is displayed in GUI
    d. name is the interface option as displayed in nipype documentation
- Naming conventions
  1. Nodes starting with a capital F (i.e. Fsmooth) are user defined functions with the ability for dynamic function inputs, and letters after F must be lowercase
    a. Function parameters must follow this naming convention:
      i. For variables only affecting the function and not internal nodes -> all lower case, one word (no underscore)
      ii. For nodes inside the function -> nodename_nipypeinterfaceparametername (similar to the above section, nipypeinterfaceparametername is copied exactly from nipype including underscores)
    b. To enable dynamic parameter assignment, there must be a linebreak before the workflow is run inside a function (i.e. randomline \n\n workflow.run())
  2. Entries in json file starting with ~construct~ (i.e. ~construct~Finfo) indicate a dictionary will be constructed from the provided information, and passed to the function (i.e. Finfo) in the format of dictionaryname_parametername
  3. Entries starting with ! (i.e. !correction) indicate that the value will be copied from another parameter as defined in default_links.json
    a. node_to_add and node_to_copy indicate the name of the parameter to add, and the node the values will be copied from, respectively
      i. An additional entry of on_off can be used to modify the above, where the parameter being copied is dependent on the node specified with on_off (i.e. true -> use value, false -> no action)
    b. verify and values together alter the construction of a dictionary, with values specifying the scenario in which verify is added to the dictionary
    c. node_to_edit is used if there are mutually exclusive options that must be handled
      i. on_off is the node which controls mutual exclusivity, switch indicates the value of on_off that is mutually exclusive with node_to_edit
  4. Naming conventions MUST be followed as they are what enable the dynamic creation of custom pipelines which share data as long as possible, and don't force an analysis of all permutations of multiverse options

# Resource Allocation
 - On compute canada ~4 days for 50 subjects x 200 pipelines with 200 CPUs, 6gb RAM per CPU
 - Generates a lot of data, peaking at ~0.83GB per subject per pipeline
   a. Running with debug set to false will delete files once they are no longer needed by the workflow
     i. Saves a LOT of space, but if the analysis fails, it cannot be rerun from where it failed, and will restart from the beginning (i.e. progress is lost)
     ii. As a consequence, give a buffer when requesting run time, as computation time will be wasted if program fails prior to exiting

  
