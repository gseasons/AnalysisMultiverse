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

processed = '/scratch/processed'
    
if len(sys.argv) > 3:
    batch = sys.argv[1]
    task = sys.argv[2]
    profile = sys.argv[3]
    
save_dir = processed + '/checkpoints_' + task + '_batch_' + str(batch)
config.set("execution", "crashdump_dir", save_dir)
    
working_dir = '/scratch/{task}_working_dir_{batch}'.format(task=task, batch=batch)

with open(processed + '/reproducibility/' + task + '_workflow_' + batch + '.pkl', 'rb') as wf:
    workflow = pickle.load(wf)
    
workflow.run(plugin='IPython', plugin_args={'profile': profile, 'task': task, 'batch': batch})

with open('/code/multiverse/configuration/general_configuration.pkl', 'rb') as f:
    configuration = pickle.load(f)

if os.path.exists(save_dir) and not glob.glob(save_dir + '/crash-*'):
    os.rename(save_dir, save_dir + '_done')

if os.path.exists(working_dir) and not configuration['debug'] and not glob.glob(save_dir + '/crash-*'):
    shutil.rmtree(working_dir)