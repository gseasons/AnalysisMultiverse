#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 12:28:44 2021

@author: grahamseasons
"""
from analysis_pipeline import *
from os.path import join as opj
from bids.layout import BIDSLayout
from os.path import join as opj
import os

def main():
    data_dir = "/Users/grahamseasons/fMRI/test_data/ds000114/"
    exp_dir = '/Volumes/NewVolume/ds000114-mult' #right'#'allsubs'
    working_dir = 'working_dir'
    
    layout = BIDSLayout(data_dir)
    subjects = layout.get_subjects()
    subjects = ['02','03']#['01','02','03','10']
    #subjects = ['02','03','04','05','07','08','09']
    types = layout.get_datatypes()
    sessions = layout.get_sessions()
    runs = layout.get_runs()
    tasks = layout.get_tasks()
    tasks = ['fingerfootlips']
    
    #GA SETUP
    num_pipelines = 1
    
    if 'anat' and 'func' in types:
        go = analysis(exp_dir, working_dir, data_dir).construct(subjects, tasks, sessions, runs, num_pipelines)
        
        #GA LOOP -> these will all be in dictionaries and assigned as outputs from a function
        go.inputs.inputnode.pipeline_names = ['test']#, 'test2']
        go.inputs.inputnode.frac_mask = [0.3]#, 0.3]
        go.inputs.inputnode.discard = [4]#, 4]
        go.inputs.inputnode.dof_mc = [6]#, 6]#not connected
        go.inputs.inputnode.fwhm = [4]#, 8]
        go.inputs.inputnode.cost_mc = ['normcorr']#, 'normcorr']#not connected
        go.inputs.inputnode.bet_frac = [0.5]#, 0.5]
        go.inputs.inputnode.robust = [True]#, True]
        go.inputs.inputnode.wm_thresh = [0.5]#, 0.5]
        go.inputs.inputnode.dof_f = [6]#, 6]
        go.inputs.inputnode.bbr_type = ['signed']#, 'signed']#temporarily removed
        go.inputs.inputnode.interp = ['spline']#, 'spline']
        go.inputs.inputnode.cost = ['mutualinfo']#, 'mutualinfo']#temporarily removed
        go.inputs.inputnode.bins = [640]#, 640]#temporarily removed
        go.inputs.inputnode.iso = [4]#, 4]
        go.inputs.inputnode.bbr = [True]#, True]
        go.inputs.inputnode.susan = [True]#, True]
        go.inputs.inputnode.HP = [128]#, 140]
        go.inputs.inputnode.film_thresh = [1000]#, 1000]#not connected
        go.inputs.inputnode.serial_cor = [True]#, True]
        go.inputs.inputnode.base_switch = [False]#, False]
        go.inputs.inputnode.gamma = [False]#, False]
        go.inputs.inputnode.dict_opt = ['gammasigma']#, 'gammasigma']#connected but unused
        go.inputs.inputnode.base_val = [3]#, 3]#connected but unused
        go.inputs.inputnode.mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
        go.inputs.inputnode.mode_l2 = ['fe']#, 'fe']
        go.inputs.inputnode.mode_l3 = ['flame1']#, 'flame1']
        go.inputs.inputnode.method = ['fwe']#, 'cluster']#'fdr']#'cluster']
        go.inputs.inputnode.p = [0.025]#, 0.025]
        go.inputs.inputnode.connectivity = [26]#, 26]
        go.inputs.inputnode.z_thresh = [3.2]#, 2.6]
        go.run()
        #grab output pipeline scores, feed into GA class
        
        #done desired number of pipelines, data output
        A = 3
        
if __name__ == "__main__":
    main()