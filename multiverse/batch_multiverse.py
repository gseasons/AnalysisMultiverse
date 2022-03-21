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
import getpass
import subprocess
from functions import organize

processed = '/scratch/processed'
    
if len(sys.argv) > 3:
    batch = sys.argv[1]
    task = sys.argv[2]
    profile = sys.argv[3]
    
working_dir = '/scratch/{task}_working_dir_{batch}'.format(task=task, batch=batch)

with open(processed + '/reproducibility/' + task + '_workflow_' + batch + '.pkl', 'rb') as wf:
    workflow = pickle.load(wf)
    
workflow.run(plugin='IPython', plugin_args={'profile': profile, 'task': task, 'batch': batch})

if os.path.exists(working_dir):
    shutil.rmtree(working_dir)

user = getpass.getuser()
sq = str(subprocess.check_output(['squeue', '--user={0}'.format(user)]))
out_frame = processed + '/{task}.pkl'.format(task=task)

if sq.count('batch.sh') == 1:
    paths = organize(task, out_frame)
    
    #INSERT DATA PROCESSING HERE - TO BE WRITTEN ONCE WE HAVE TEST DATA
    
