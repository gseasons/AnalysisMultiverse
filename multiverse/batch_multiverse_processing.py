#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:12:07 2022

@author: grahamseasons
"""
import sys
from functions import organize

processed = '/scratch/processed'
    
if len(sys.argv) > 2:
    task = sys.argv[1]
    sq = sys.argv[2]
    
out_frame = processed + '/{task}.pkl'.format(task=task)

if sq.count('batch.sh') == 1:
    paths = organize(task, out_frame)
    
    #INSERT DATA PROCESSING HERE - TO BE WRITTEN ONCE WE HAVE TEST DATA