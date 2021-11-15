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
from collections import Counter

def define_paths(container, dictionary, indexes):
    #keys, values = zip(*dictionary.items())
    #will probably change this into comparing hashes of dictionaries or something: 
        #https://www.doc.ic.ac.uk/~nuric/coding/how-to-hash-a-dictionary-in-python.html
    out_dic = {}
    link = {}
    old_x = len(indexes)
    if old_x >= 1:
        old_x -= 1
    change = len(indexes) - 1
    for i, vals in enumerate(indexes):
        if isinstance(vals, np.ndarray):
            out_dic[i] = {'id': vals}
        else:
            out_dic[i] = {'id': np.array(range(vals))}
    out_dic['const'] = {}
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
                            
                        vals, ind = np.unique(container, return_inverse=True, axis=1)
                        
                        if True:#vals.shape != container.shape:
                            index = [np.where((vals[:,i].reshape(-1,1) == container).sum(axis=0) == container.shape[0])[0] for i in range(vals.shape[1])]
                            ind = index.copy()
                            for arr in index:
                                ind[arr.min()] = arr
                            index_ = ind
                            gen = np.unique(tup[1])
                            for k in out_dic:
                                if type(k) != int:#'const':
                                    break
                                if len(gen) == 1:
                                    if key not in out_dic['const']:
                                        out_dic['const'][key] = {}
                                    if subkey not in out_dic['const'][key]:
                                        out_dic['const'][key][subkey] = {}
                                        
                                    out_dic['const'][key][subkey][tup[0]] = tup[1][0]
                                    continue
                                    
                                for j, x in enumerate(index_):
                                    if 'id' in out_dic[k]:
                                        if len(np.intersect1d(x, out_dic[k]['id'])) != len(x):
                                            continue
                                    if j not in out_dic[k]:
                                        out_dic[k][j] = {}
                                    if key not in out_dic[k][j]:
                                        out_dic[k][j][key] = {}
                                    if subkey not in out_dic[k][j][key]:
                                        out_dic[k][j][key][subkey] = {}
                                    if 'id' not in out_dic[k][j]:
                                        out_dic[k][j]['id'] = [-1]
                                        
                                    if not np.array_equiv(out_dic[k][j]['id'], x):
                                        if not np.array_equiv(out_dic[k][j]['id'], [-1]):
                                            cx = Counter(out_dic[k][j]['id'])
                                            cid = Counter(x)
                                            out = min(sorted((cx - cid).elements()))
                                            link = {out: {j: [key, subkey]}}
                                        elif 'id' in out_dic[k]:
                                            if not np.array_equiv(out_dic[k][j]['id'], [-1]) or (min(x) == min(out_dic[k]['id']) and len(x) < len(out_dic[k]['id'])):
                                                #cx = Counter(out_dic[k]['id'])
                                                #cid = Counter(x)
                                                #out = min(sorted((cx - cid).elements()))
                                                for out in out_dic[k]['id']:
                                                    link[out] = {j: [key, subkey]}
                                                #link = {out: {j: [key, subkey]}}
                                    
                                    if j > old_x:
                                        out_dic[k][j]['link'] = link[j]
                                    
                                    out_dic[k][j]['id'] = x
                                    out_dic[k][j][key][subkey][tup[0]] = tup[1][j]
                                if k == change:#list(out_dic.keys())[-2]:
                                    old_x = j
                                
                            link = {}
                            out_dic[subkey + '_' + tup[0]] = tup[1]
                    
    un, ind = np.unique(container, return_inverse=True, axis=1)
    index = [np.where((un[:,i].reshape(-1,1) == container).sum(axis=0) == container.shape[0])[0] for i in range(un.shape[1])]
    for arr in index:
        ind[arr] = arr.min()
        
    return container, list(ind), out_dic, index_
            

def main():
    data_dir = '/Volumes/NewVolume/super_agers' #"/Users/grahamseasons/fMRI/test_data/ds000114/"
    exp_dir = '/Volumes/NewVolume/sup_test' #'/Volumes/NewVolume/ds000114-datasink' #right'#'allsubs'
    working_dir = 'working_dir'
    out_dir = exp_dir + '/processed' #derivatives probably output to same directory as data
    
    layout = BIDSLayout(data_dir)
    subjects = ['02','03']#['01','02','03','10']
    #subjects = ['02','03','04','05','07','08','09']
    tasks = layout.get_tasks()
    #tasks = ['fingerfootlips']
    
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
    resting = [{'atlas': 'harvard', 'seed': 'cing post', 'seed_thr': 25, 'x': 45, 'y': 74, 'z': 51, 'radius': 5,
                'CSF': True, 'GLOBAL': True, 'CSF_thresh': 0.5}, 
               {'atlas': 'ROI', 'seed': '', 'seed_thr': 0, 'x': 45, 'y': 74, 'z': 51, 'radius': 5, 
                'CSF': True, 'GLOBAL': True, 'CSF_thresh': 0.5}]
    
    group_num = 10
    
    for task in tasks:
        subjects = layout.get_subjects(task=task)[0:10]
        types = layout.get_datatypes()
        sessions = layout.get_sessions(task=task)
        runs = layout.get_runs(task=task)
        
        #subjects = ['02','03']#['01','02','03','10']
        #subjects = ['02','03','04','05','07','08','09']
        
        preproc = {}
        spatial_norm = {}
        level1 = {}
        level2 = {}
        level3 = {}
        correction = {}
        container = np.zeros([num_pipelines+2])
        #IN FINAL PRODUCT PREPROC IS PASSED IN, BUT EVERY OTHER INPUT IS SET WITH INPUTNODE CONNECTIONS FROM HERE
        if 'anat' and 'func' in types:                                                        
            preproc['improcess'] = {'ExtractROI': [('discard', discard+[4]*2)],
                                    #'SliceTimer': [],
                                    #'MCFLIRT': [],
                                    'smoother': [('susan', susan+[True]*2), ('fwhm', fwhm+[4]*2), ('frac_mask', frac_mask+[0.5]*2)],
                                    'PIPELINE': [('pipeline_names', pipeline_names)],
                                    }
            egg = {}
            egg['coreg'] = {'BET': [('bet_frac', bet_frac + [0.2, 0.3]), ('robust', robust+[False, False])],
                                #'FAST': [],
                                'registration': [('bbr', bbr), ('wm_thresh', wm_thresh), ('dof_f', dof_f), ('bbr_type', bbr_type), ('interp', interp),
                                                 ('iso', iso), ('cost', cost), ('bins', bins), ('warplater', warp_post_feat)],
                                'PIPELINE': [('pipeline_names', pipeline_names)],
                                }
            
            container, ind, out_dic, index_ = define_paths(container, preproc, [num_pipelines+2])
            container, ind, out_dic, index_ = define_paths(container, egg, index_)
            
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
            
            go = analysis(exp_dir, working_dir, data_dir, out_dir).construct(subjects, sessions, runs, task, group_num, pipeline_ind, preproc, spatial_norm, level1, level2, level3, correction)
            go.inputs.inputnode.mask = mask
            go.inputs.inputnode.task = task
            go.inputs.inputnode.max_groups = group_num                                  
            
            #GA LOOP -> these will all be in dictionaries and assigned as outputs from a function
            go.run(plugin='MultiProc')
        #grab output pipeline scores, feed into GA class
        
        #done desired number of pipelines, data output
        A = 3
        
if __name__ == "__main__":
    main()