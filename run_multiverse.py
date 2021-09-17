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
import numpy as np

def define_paths(container, dictionary):
    #keys, values = zip(*dictionary.items())
    #will probably change this into comparing hashes of dictionaries or something: 
        #https://www.doc.ic.ac.uk/~nuric/coding/how-to-hash-a-dictionary-in-python.html
    for i, key in enumerate(dictionary):
        if type(dictionary[key]) == list and key != 'PIPELINE':
            for tup in dictionary[key]:
                if type(tup[1][0]) != dict:
                    container = np.vstack((container, np.array(tup[1])))
                else:
                    to_stack = np.array([list(subdict.values()) for subdict in tup[1]]).transpose()
                    container = np.vstack((container, to_stack))
        elif type(dictionary[key]) == dict:
            for subkey in dictionary[key]:
                if type(dictionary[key][subkey]) == list and subkey != 'PIPELINE':
                    for tup in dictionary[key][subkey]:
                        if type(tup[1]) == dict:
                            for subsubkey in tup[1]:
                                container = np.vstack((container, np.array(tup[1][subsubkey])))
                        else:
                            container = np.vstack((container, np.array(tup[1])))
        
    un, ind = np.unique(container, return_inverse=True, axis=1)
    index = [np.where((un[:,i].reshape(-1,1) == container).sum(axis=0) == container.shape[0])[0] for i in range(un.shape[1])]
    for arr in index:
        ind[arr] = arr.min()
        
    return container, list(ind)
            

def main():
    data_dir = "/Users/grahamseasons/fMRI/test_data/ds000114/"
    exp_dir = '/Volumes/NewVolume/ds000114-iter' #right'#'allsubs'
    working_dir = 'working_dir'
    
    layout = BIDSLayout(data_dir)
    subjects = ['02','03']#['01','02','03','10']
    #subjects = ['02','03','04','05','07','08','09']
    tasks = layout.get_tasks()
    tasks = ['fingerfootlips']
    
    #GA SETUP
    num_pipelines = 2
    #GA DEFINE PIPELINES -> FINAL VERSION WILL OUTPUT DICTIONARIES
    pipeline_names = ['test', 'test2'] #-> WILL BE USED FOR WRITING TO DATASINK, NEEDS TO BE PASSED TO EVERY PLACE AN OUTPUT IS KEPT/SAVED
    frac_mask = [0.3, 0.3]
    discard = [4, 4]
    dof_mc = [6, 6]#not connected
    fwhm = [4, 8]
    cost_mc = ['normcorr', 'normcorr']#not connected
    bet_frac = [0.5, 0.5]
    robust = [True, True]
    wm_thresh = [0.5, 0.5]
    dof_f = [6, 6]
    bbr_type = ['signed', 'signed']#temporarily removed
    interp = ['spline', 'spline']
    cost = ['mutualinfo', 'mutualinfo']#temporarily removed
    bins = [640, 640]#temporarily removed
    iso = [4, 4]
    bbr = [True, True]
    susan = [True, True]
    HP = [128, 140]
    film_thresh = [1000, 1000]#not connected
    serial_cor = [True, True]
    base_switch = [False, False]
    gamma = [False, False]
    dict_opt = ['gammasigma', 'gammasigma']#connected but unused
    base_val = [3, 3]#connected but unused
    mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
    mode_l2 = ['fe', 'fe']
    mode_l3 = ['flame1', 'flame1']
    method = ['fwe', 'cluster']#'fdr']#'cluster']
    p = [0.025, 0.025]
    connectivity = [26, 26]
    z_thresh = [3.2, 2.6]
    warp_post_feat = [True, True]
    resting = [{'atlas': 'harvard', 'seed': 'cing post', 'seed_thr': 25}]*2
    
    for task in tasks:
        subjects = layout.get_subjects(task=task)
        types = layout.get_datatypes(task=task)
        sessions = layout.get_sessions(task=task)
        runs = layout.get_runs(task=task)
        
        subjects = ['02','03']#['01','02','03','10']
        #subjects = ['02','03','04','05','07','08','09']
        
        preproc = {}
        spatial_norm = {}
        level1 = {}
        level2 = {}
        level3 = {}
        correction = {}
        container = np.zeros([num_pipelines])
        #IN FINAL PRODUCT PREPROC IS PASSED IN, BUT EVERY OTHER INPUT IS SET WITH INPUTNODE CONNECTIONS FROM HERE
        if 'anat' and 'func' in types:                                                        
            preproc['improcess'] = {'ExtractROI': [('discard', discard)],
                                    #'SliceTimer': [],
                                    #'MCFLIRT': [],
                                    'smoother': [('susan', susan), ('fwhm', fwhm), ('frac_mask', frac_mask)],
                                    'PIPELINE': [('pipeline_names', pipeline_names)],
                                    }
            
            preproc['coreg'] = {'BET': [('bet_frac', bet_frac), ('robust', robust)],
                                #'FAST': [],
                                'registration': [('bbr', bbr), ('wm_thresh', wm_thresh), ('dof_f', dof_f), ('bbr_type', bbr_type), ('interp', interp),
                                                 ('iso', iso), ('cost', cost), ('bins', bins), ('warplater', warp_post_feat)],
                                'PIPELINE': [('pipeline_names', pipeline_names)],
                                }
            
            container, ind = define_paths(container, preproc)
            
            spatial_norm['get'] = {}
            spatial_norm['apply'] = {'WARP': [('needwarp', warp_post_feat)],
                                     'PIPELINE': [('pipeline_names', pipeline_names)],
                                     }
            
            container, ind = define_paths(container, spatial_norm['apply'])
            spatial_norm['apply']['ind'] = [('ind', ind)]
            
            
            #FINAL GOAL -> USE DICTIONARIES AS ITERABLES i.e. EACH WILL BE A LIST OF DICTIONARIES, SO THAT DICTIONARIES CAN BE PASSED TO FUNCTIONS, AND FED INTO NODES DIRECTLY
            
            level1 = {'session_info': [('HP', HP), ('warp_post_feat', warp_post_feat), ('resting', resting)],
                      'l1d': [('gamma', gamma), ('dict_opt', dict_opt), ('base_val', base_val), ('base_switch', base_switch), ('serial_cor', serial_cor),
                              ('discard', discard)],
                      'PIPELINE': [('pipeline_names', pipeline_names)],
                      }
            
            container, ind = define_paths(container, level1)
            level1['ind'] = [('ind', ind)]
            
            #TODO: DELETE ALL CONNECTIONS AND INPUTNODE MENTIONS OF SUBCATEGORIES
            
            level2 = {'FLAMEO': [('run_mode', mode_l2)],
                      'PIPELINE': [('pipeline_names', pipeline_names)],
                      }
            
            container, ind = define_paths(container, level2)
            level2['ind'] = [('ind', ind)]
            pipeline_ind = [('ind', ind)]
            
            level3 = {'FLAMEO': [('run_mode', mode_l3)],
                      'PIPELINE': [('pipeline_names', pipeline_names)],
                      }
            
            container, ind = define_paths(container, level3)
            level3['ind'] = [('ind', ind)]
            
            correction['decision'] = {'DEC': [('method', method), ('p', p), ('connectivity', connectivity), ('z_thresh', z_thresh), ('pipeline_names', pipeline_names)],
                                      'PIPELINE': [('pipeline_names', pipeline_names)],
                                      }
            
            container, ind = define_paths(container, correction)
            correction['decision']['ind'] = [('ind', ind)]
            
            go = analysis(exp_dir, working_dir, data_dir).construct(subjects, sessions, runs, pipeline_ind, preproc, spatial_norm, level1, level2, level3, correction)
            go.inputs.inputnode.mask = mask
            go.inputs.inputnode.task = task                                          
            
            #GA LOOP -> these will all be in dictionaries and assigned as outputs from a function
            go.run()
        #grab output pipeline scores, feed into GA class
        
        #done desired number of pipelines, data output
        A = 3
        
if __name__ == "__main__":
    main()