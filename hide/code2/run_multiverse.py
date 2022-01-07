#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:41:49 2021

@author: grahamseasons
"""
import pygad as pg
import pickle
import numpy as np
import re, os, random
from os.path import join as opj
#from updated.functions import define_paths
from functions import define_paths, generate_dictionaries
#from updated.analysis_pipeline import analysis
from analysis_pipeline import analysis
from bids.layout import BIDSLayout
from nipype.interfaces.base import Undefined
#from updated.functions import generate_dictionaries
#from functions import generate_dictionaries
#from nipype.interfaces.traits import _Undefined
from nipype import config
#from nipype.utils.profiler import log_nodes_cb
#config.enable_debug_mode()
config.set("execution", "hash_method", "content")
config.set("execution", "remove_node_directories", "true")
#config.enable_resource_monitor()
exp_dir = '/scratch'#'/Volumes/NewVolume/sup_pre_full'
working_dir = 'working_dir'
data_dir = '/data'#'/Volumes/NewVolume/super_agers'
out_dir = exp_dir + '/processed'
mask = opj(os.getenv('FSLDIR'), 'data/standard/MNI152_T1_2mm_brain.nii.gz')

links = {'preprocess': {
         'decision_mc_mean': 'mcflirt_mean_vol',
         'Fregress_glm_des_norm': 'Fregress_glm_dat_norm',
         #'Fregistration_wmthresh': ['Fregistration_wmthresh', 'Fregistration_bbr'],
         'Fregistration_regbbr_interp': ['Fregistration_regpre_interp', 'Fregistration_bbr'],
         'Fregistration_regbbr_no_resample': ['Fregistration_regpre_no_resample', 'Fregistration_bbr'],
         'Fregistration_applywarp_interp': 'Fregistration_regpre_interp',
         'Fregistration_applywarp_no_resample': 'Fregistration_regpre_no_resample',# 'Fregistration_bbr'],
         'invwarp_warplater': 'Fmni_warplater',
         'invwarp_concatenate': 'Fregistration_concatenate',
         'Fmni_concatenate': 'Fregistration_concatenate',
         },
         'level1': {
         'seedinfo': ['atlas', 'data'],
         'atlas': ['atlas', 'data'],
         'coords': ['ROI'],
         'radius': ['ROI'],
         'k': ['data'],
         'kcc': ['data'],
         'lp': ['data'],
         'hp': ['data'],
         'Finfo_warppostfeat': 'Fmni_warplater',
         'Finfo_concatenate': 'Fregistration_concatenate',
         #'Finfo_realignregress': 'Fregress_realignregress',
         'correction_discard': 'extract_t_min',
         'applywarpcopes_interp': 'Fmni_warped_interp',
         'applywarpvarcopes_interp': 'Fmni_warped_interp',
         'applywarpbold_interp': 'Fmni_warped_interp',
         'ret_needwarp': 'Fmni_warplater',
         'ident_needwarp': 'Fmni_warplater',
         'correction_discard': 'extract_t_min',
         },
         'level2': {
         'flameo_infer_outliers': ['flameo_infer_outliers', 'flameo_run_mode'],
         },
         'level3': {
         'flameo_infer_outliers': ['flameo_infer_outliers', 'flameo_run_mode', 'fe'],#AVOID LAST ONE
         },
         'correction': {
         'pthreshold': ['clust'],
         'zthreshold': ['clust'],
         'connectivity': ['clust'],
         'p': ['fdr', 'fwe'],
         },
         }
#NOTE: #** means input cannot be removed
genes = [{'dilateref_kernel_size': [0]},#[0, 4, 8]},
         {'dilateref_kernel_shape': [5],
              0: '3D',
              1: '2D',
              2: 'box',
              3: 'boxv',
              4: 'gauss',
              5: 'sphere'},
         {'extract_t_min': [10]},
         {'mcflirt_interpolation': [0],#, 1, 2],
             0: 'spline',
             1: 'nn',
             2: 'sinc'},
         {'mcflirt_mean_vol': [0],#, 1],
              0: False,
              1: True},
         {'decision_slice_correct': [0],#, 1],
              0: False,
              1: True},
         #
         {'prelim_interp': [0],#, 1, 2, 3],
              0: 'trilinear',
              1: 'nearestneighbour',
              2: 'sinc',
              3: 'spline'},
         {'prelim_no_resample': [0],
              0: False,
              1: True},
         {'dilatebrain_kernel_size': [0]},#, 4, 8]},
         {'dilatebrain_kernel_shape': [5],
              0: '3D',
              1: '2D',
              2: 'box',
              3: 'boxv',
              4: 'gauss',
              5: 'sphere'},
         #DILATEREF USED TO BE HERE
         {'warp_warp_resolution': [0],#, 1],
              0: (10,10,10),
              1: (20,20,20)},
         #EVERYTHING ABOVE BET USED TO BE HERE
         {'Fregistration_bbr': [1],#**
              0: False,
              1: True},
         {'Fregistration_regpre_interp': [0],#, 1, 2, 3],
              0: 'trilinear',
              1: 'nearestneighbour',
              2: 'sinc',
              3: 'spline'},
         {'Fregistration_regpre_no_resample': [0], #Fregistration_regpre_no_resample
              0: False,
              1: True},
         {'Fregistration_applywarp_': [1], #NOTE: if setting a node input inside function to Undefined -> use empty string
              0: {'apply_isoxfm': 4, 'apply_xfm': '', 'uses_qform': ''},
              1: {'apply_isoxfm': '', 'apply_xfm': True, 'uses_qform': ''}},
         {'Fregistration_concatenate': [1],
              0: False,
              1: True},
         #{'Fregistration_wmthresh': {'low': 0.2, 'high': 0.8, 'step': 0.1}},#**
         {'Fmni_warplater': [1],
              0: False,
              1: True},
         {'Fmni_warped_interp': [1],#, 1],#, 2, 3],#**
              0: 'nn',
              1: 'trilinear',
              2: 'sinc',
              3: 'spline'},
         {'!invwarp_': [0],
              0: 'LINK'},
         {'Fregress_CSF': [1],#, 1],
              0: False,
              1: True},
         {'Fregress_WM': [1],#, 1],
              0: False,
              1: True},
         {'Fregress_GLOBAL': [1],#, 1],
              0: False,
              1: True},
         {'Fregress_glm_dat_norm': [1],
              0: False,
              1: True},
         {'Fregress_glm_demean': [1],
              0: False,
              1: True},
         {'Fregress_realignregress': [1],#, 1],
              0: False,
              1: True},
         {'Fsmooth_susan': [1],#**
              0: False,
              1: True},
         {'Fsmooth_fwhm': [6]}, #range(2, 8, 2)},#13, 10)},#2)},#**
         #ART
         {'end_preprocess': 'level1'},
         #Finfo_rest OPTIONS WILL NEED TO BE PARTIALLY USER DEFINED -> gui to ask?
         {'~construct~Finfo_rest_type': [0, 1, 2],
              0: 'atlas',
              1: 'ROI',
              2: 'data'},
         {'~construct~Finfo_rest_seedinfo': {'low': 0.6, 'high': 1.4},#{'low': -0.5, 'high': 2.5},
              0: [('cing post', 1)], #can include atlas name with seed info -> user specified
              1: [('cing post', 0, 1)], #different variations on same mask, or list masks targeting different networks
              2: [('cing post', 1)]},
         {'~construct~Finfo_rest_coords': [0],
              0: [(44, 35, 44)],#[(),(),()]
              1: [(44, 37, 49)],
              2: [(47, 38, 56)]},
         {'~construct~Finfo_rest_radius': [6]},#range(3, 11)},
         {'~construct~Finfo_rest_atlas': [0],
              0: 'harvard'},
         {'~construct~Finfo_rest_k': [2],#0, 1, 2],
              0: 'faces',
              1: 'edges',
              2: 'vertices'},
         {'~construct~Finfo_rest_kcc': [0.5]},#, 0.6, 0.7]},
         {'~construct~Finfo_rest_lp': [0.01]},
         {'~construct~Finfo_rest_hp': [0.1]},
         {'Finfo_HP': [128]},#range(50, 150)},
         {'!correction_': [0],
              0: 'LINK'},
         {'l1d_bases': [3],
              0: {'gamma': {'derivs': True}},
              1: {'gamma': {'derivs': False}},
              2: {'dgamma': {'derivs': True}},
              3: {'dgamma': {'derivs': False}},
              4: {'custom': {'bfcustompath': '${FSLDIR}/etc/default_flobs.flobs'}}},
         {'l1d_model_serial_correlations': [0, 1],
              0: False,
              1: True},
         {'!applywarpcopes_': [0],
              0: 'LINK'},
         {'!applywarpvarcopes_': [0],
              0: 'LINK'},
         {'!applywarpbold_': [0],
              0: 'LINK'},
         {'!ret_': [0],
              0: 'LINK'},
         {'!ident_': [0],
              0: 'LINK'},
         {'end_level1': 'level2'},
         {'flameo_run_mode': [2],#0, 1, 2, 3],
              0: 'fe',
              1: 'ols',
              2: 'flame1',
              3: 'flame12'},
         {'flameo_infer_outliers': [0],#, 1],
              0: False,
              1: True},
         {'end_level2': 'level3'},
         {'flameo_run_mode': [0, 1, 2, 3],
              0: 'fe',
              1: 'ols',
              2: 'flame1',
              3: 'flame12'},
         {'flameo_infer_outliers': [0, 1],
              0: False,
              1: True},
         {'end_level3': 'correction'},
         {'~construct~correct_cor_method': [0, 1, 2],
              0: 'fdr',
              1: 'fwe',
              2: 'clust'},
         {'~construct~correct_cor_connectivity': [26]},
         {'~construct~correct_cor_zthreshold': [2, 2.3, 2.7, 3.1, 3.5, 4]},
         {'~construct~correct_cor_pthreshold': [0.05]},
         {'~construct~correct_cor_p': [0.1, 0.05, 0.01, 0.005, 0.001]},
         {'end_correction': 'end'},
         ]

#THESE WILL BE BOUNDS AROUND SINGLE VALUE IF GIVEN -> CHANGE SO DOESN'T GO UNDER 5 or over 95
wiggle = 10
split_half = False

map_genes = {}
#WORKS FOR NOW
def fitness_func(solution, solution_idx):
    #GLOB FOR PROPER dist VALUE FOR SOLUTION, READ VALUE IN -> FITNESS
    #MIGHT HAVE TO DO SOMETHING TRICKY TO FIND RIGHT PIPELINE AS IDX IS ONLY FOR SPECIFIC POPULATION
    #AND DOESN'T APPEAR TO BE ACCESS TO GENERATION NUMBER -> could right folders as index number, then
    #use that to search, and rewrite filename after that
    import random
    return random.randint(0, 90)


def on_pop_gen(ga): #on_generation, on_start, check to make sure not rerunning pipelines
    gen = ga.generations_completed
    generation = gen
    pop = ga.population
    params = pop.transpose()
    pipeline = generation * pop.shape[0]
    print('\n\n\n\n')
    print(str(pipeline))
    layout = BIDSLayout(data_dir)
    tasks = layout.get_tasks()
    for task in tasks:
        subjects = layout.get_subjects(task=task)
        subjects.sort()
        subjects = subjects[0:1]
        types = layout.get_datatypes()
        sessions = layout.get_sessions(task=task)
        runs = layout.get_runs(task=task)
        
        if sessions or runs:
            multiscan = True
        else:
            multiscan = False
        #ALTER SO THAT CHECKS AGAINST DICTIONARY TO SEE WHAT VALUES ARE EXCLUSIVE FOR RESTING STATE -> IGNORES THEM IN TASK
        #NAME AND DEFAULT VALUE
        master, expand_inputs = generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle)
        #output = open('master.pkl', 'wb')
        #pickle.dump(master, output)
        #output.close()
        #output = open('expand.pkl', 'wb')
        #pickle.dump(expand_inputs, output)
        #output.close()
        #output = open('population.pkl', 'wb')
        #pickle.dump(params, output)
        #output.close()
        with open('/scratch/master.pkl', 'rb') as f:
            master = pickle.load(f)
        with open('/scratch/expand.pkl', 'rb') as f:
            expand_inputs = pickle.load(f)
        if 'anat' in types and 'func' in types:
            pipelines = analysis(exp_dir, working_dir, data_dir, out_dir)
            pipelines = pipelines.construct(subjects, sessions, runs, task, pipeline, master, expand_inputs, split_half)
            pipelines.inputs.inputnode.mask = mask
            pipelines.inputs.inputnode.task = task
            pipelines.run(plugin='MultiProc')#, plugin_args={'status_callback': log_nodes_cb})
    
    A = 3
    import sys
    sys.exit()
                            
                    

def main():
    num_generations = 5
    num_parents_mating = 4
    
    gene_space = []
    dummy = 0
    num_genes = len(genes)
    
    for i, gene in enumerate(genes):
        map_genes[i] = gene
        if 'end' in list(gene.keys())[0]:
            dummy += 1
            continue
        
        vals = list(gene.values())[0]
        if type(vals) == list:
            vals = [int(val) if type(val) == bool else val for val in vals]
            
        gene_space.append(vals)
    
    num_genes -= dummy
    
    ga = pg.GA(num_generations=num_generations,
               num_parents_mating=num_parents_mating,
               fitness_func=fitness_func,
               on_start=on_pop_gen,
               on_generation=on_pop_gen,
               sol_per_pop=8,
               
               gene_type=[float, 3],
               
               num_genes=num_genes,
               gene_space=gene_space,
               parent_selection_type='sss',
               keep_parents=-1,
               crossover_type="single_point",
               mutation_type="random",
               mutation_probability=0.2,
               )
    ga.run()
    
if __name__ == "__main__":
    main()

