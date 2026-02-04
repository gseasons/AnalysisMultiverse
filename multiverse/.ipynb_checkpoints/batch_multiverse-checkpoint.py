#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 13:58:49 2022

@author: grahamseasons
"""
import pickle
import sys
import shutil
import os
import glob
from nipype import config
import ipyparallel #Mah
import subprocess  #Mah

print("we are in the batch_multiverse.py | MAH") #mah

processed = '/scratch_dir/processed'
    
if len(sys.argv) > 3:
    print('We are in the if len(sys.argv) > 3 |Mah') #Mah 
    batch = sys.argv[1]
    task = sys.argv[2]
    profile = sys.argv[3]

# print('MAH| batch : ', sys.argv[1], ' task: ',sys.argv[2] , ' profile: ', sys.argv[3]) #Mah
    
save_dir = processed + '/reproducibility/checkpoints_' + task + '_batch_' + str(batch)
if not os.path.exists(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    
config.set("execution", "crashdump_dir", save_dir)

with open('/code/multiverse/configuration/general_configuration.pkl', 'rb') as f:
    configuration = pickle.load(f)

if not configuration['debug']:
    config.set("execution", "remove_node_directories", "true")
    
working_dir = '/scratch_dir/{task}_working_dir_{batch}'.format(task=task, batch=batch)

with open(processed + '/reproducibility/' + task + '_workflow_' + batch + '.pkl', 'rb') as wf:
    workflow = pickle.load(wf)
    # print('workflow.client_args: ', workflow.client_args) #Mah


#----------------------------------
import os
print("MAH| Checking if IPython profile directory exists:",
      os.path.exists(os.path.expanduser(f"~/.ipython/profile_{profile}")))

# print('MAH| print :os.getcwd(): ', os.getcwd())

# # Run 'ls' and capture output #Mah
# result = subprocess.run(["ls", "-a"], capture_output=True, text=True)
# print('MAH| printing directories in ls:', result.stdout) #This correctly shows .ipython in this directory


# result0 = subprocess.run(["pwd"], capture_output=True, text=True)
# print('MAH| printing directories in pwd:', result0.stdout) #This correctly shows .ipython in this directory

# # Run 'ls .ipython' and capture output #Mah #but when i use this code, it says ~/.ipython doesnt exist
# home_ipython_dir = os.path.expanduser("~/.ipython")

# # Check if directory exists first
# if os.path.isdir(home_ipython_dir):
#     result = subprocess.run(["ls", home_ipython_dir], capture_output=True, text=True)
#     print("Contents of ~/.ipython:\n", result.stdout)
# else:
#     print("Directory ~/.ipython does not exist.")
#---------------------------------

workflow.run(plugin='IPython', plugin_args={'profile': profile, 'task': task, 'batch': batch})


if os.path.exists(save_dir) and not glob.glob(save_dir + '/crash-*'):
    os.rename(save_dir, save_dir + '_done')

if os.path.exists(working_dir) and not configuration['debug'] and not glob.glob(save_dir + '/crash-*'):
    shutil.rmtree(working_dir)
