#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:41:49 2021

@author: grahamseasons
"""
import pygad as pg
import numpy as np
import re, os, random
from os.path import join as opj
from updated.functions import define_paths
from updated.analysis_pipeline import analysis
from bids.layout import BIDSLayout
from nipype.interfaces.base import Undefined
#from nipype.interfaces.traits import _Undefined

exp_dir = '/Volumes/NewVolume/sup_preprocess_test_join_reg_fix'
working_dir = 'working_dir'
data_dir = '/Volumes/NewVolume/super_agers'
out_dir = exp_dir + '/processed'
mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
group_num = 10

links = {'preprocess': {
         'decision_mc_mean': 'mcflirt_mean_vol',
         'Fregress_glm_des_norm': 'Fregress_glm_dat_norm',
         #'Fregistration_wmthresh': ['Fregistration_wmthresh', 'Fregistration_bbr'],
         'Fregistration_regbbr_interp': ['Fregistration_regpre_interp', 'Fregistration_bbr'],
         'Fregistration_regbbr_no_resample': ['Fregistration_regpre_no_resample', 'Fregistration_bbr'],
         'Fregistration_applywarp_interp': 'Fregistration_regpre_interp',
         'Fregistration_applywarp_no_resample': 'Fregistration_regpre_no_resample'
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
         'invwarp_warppostfeat': 'Fmni_warplater',
         'Finfo_warppostfeat': 'Fmni_warplater',
         'correction_discard': 'extract_t_min',
         'applywarpcopes_interp': 'Fmni_warped_interp',
         'applywarpvarcopes_interp': 'Fmni_warped_interp',
         'applywarpbold_interp': 'Fmni_warped_interp',
         'ret_needwarp': 'Fmni_warplater',
         'ident_needwarp': 'Fmni_warplater'
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
         #
         #RUN SEQUENTIALLY
# =============================================================================
#          {'bet_frac': [0.5]}, #{'low': 0.15, 'high': 0.7, 'step': 0.05}},
#          {'bet_vertical_gradient': [0.2]},#{'low': -0.5, 'high': 0.5, 'step': 0.05}},
#          {'bet_': [1],#[0, 1, 2],
#               0: {'robust': True, 'reduce_bias': False, 'remove_eyes': False},
#               1: {'robust': False, 'reduce_bias': True, 'remove_eyes': False},
#               2: {'robust': False, 'reduce_bias': False, 'remove_eyes': True}},
# =============================================================================
         {'extract_t_min': [0]},
         {'mcflirt_interpolation': [0, 1, 2],
             0: 'spline',
             1: 'nn',
             2: 'sinc'},
         {'mcflirt_mean_vol': [0, 1],
              0: False,
              1: True},
         {'decision_slice_correct': [0, 1],
              0: False,
              1: True},
         {'Fregress_CSF': [0, 1],
              0: False,
              1: True},
         {'Fregress_WM': [0, 1],
              0: False,
              1: True},
         {'Fregress_GLOBAL': [0, 1],
              0: False,
              1: True},
         {'Fregress_glm_dat_norm': [1],
              0: False,
              1: True},
         {'Fregress_glm_demean': [1],
              0: False,
              1: True},
         {'Fmni_realignregress': [0, 1],
              0: False,
              1: True},
         {'Fmni_warplater': [0, 1],
              0: False,
              1: True},
         #
         {'prelim_interp': [0, 1, 2, 3],
              0: 'trilinear',
              1: 'nearestneighbour',
              2: 'sinc',
              3: 'spline'},
         {'prelim_no_resample': [0, 1],
              0: False,
              1: True},
         {'dilatebrain_kernel_size': [0, 4, 8]},
         {'dilatebrain_kernel_shape': [5],
              0: '3D',
              1: '2D',
              2: 'box',
              3: 'boxv',
              4: 'gauss',
              5: 'sphere'},
         #DILATEREF USED TO BE HERE
         {'warp_warp_resolution': [0, 1],
              0: (10,10,10),
              1: (20,20,20)},
         
         #EVERYTHING ABOVE BET USED TO BE HERE
         {'Fregistration_bbr': [1],#**
              0: False,
              1: True},
         {'Fregistration_regpre_interp': [0, 1, 2, 3],
              0: 'trilinear',
              1: 'nearestneighbour',
              2: 'sinc',
              3: 'spline'},
         {'Fregistration_regpre_no_resample': [0, 1],
              0: False,
              1: True},
         {'Fregistration_applywarp_': [1],
              0: {'apply_isoxfm': 4, 'apply_xfm': Undefined, 'uses_qform': Undefined},
              1: {'apply_isoxfm': Undefined, 'apply_xfm': True, 'uses_qform': True}},
         {'Fmni_warped_interp': [0, 1, 2, 3],#**
              0: 'nn',
              1: 'trilinear',
              2: 'sinc',
              3: 'spline'},
         #{'Fregistration_wmthresh': {'low': 0.2, 'high': 0.8, 'step': 0.1}},#**
         {'Fsmooth_susan': [1],#**
              0: False,
              1: True},
         {'Fsmooth_fwhm': range(2, 13, 2)},#**
         #ART
         {'end_preprocess': 'level1'},
         #Finfo_rest OPTIONS WILL NEED TO BE PARTIALLY USER DEFINED -> gui to ask?
         {'~construct~Finfo_rest_type': [0, 1, 2],
              0: 'atlas',
              1: 'ROI',
              2: 'data'},
         {'~construct~Finfo_rest_seedinfo': {'low': -0.5, 'high': 2.5},
              0: [('cing post', 25)], #can include atlas name with seed info -> user specified
              1: [('cing post', 5, 95)], #different variations on same mask, or list masks targeting different networks
              2: [('cing post', 50)]},
         {'~construct~Finfo_rest_coords': [0],
              0: [(50, 25, 60),(75, 30, 60),(100, 85, 62)],
              1: [(),(),()],
              2: [(),(),()]},
         {'~construct~Finfo_rest_radius': range(3, 11)},
         {'~construct~Finfo_rest_atlas': [0],
              0: 'harvard'},
         {'~construct~Finfo_rest_k': [0, 1, 2],
              0: 'faces',
              1: 'edges',
              2: 'vertices'},
         {'~construct~Finfo_rest_kcc': [0.5, 0.6, 0.7]},
         {'~construct~Finfo_rest_lp': [0.01]},
         {'~construct~Finfo_rest_hp': [0.1]},
         {'!invwarp_': [0],
              0: 'LINK'},
         {'Finfo_HP': range(50, 150)},
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
         {'end_level1': 'level2'}
         ]

#THESE WILL BE BOUNDS AROUND SINGLE VALUE IF GIVEN -> CHANGE SO DOESN'T GO UNDER 5 or over 95
edge1 = 10
edge2 = 10


map_genes = {}

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
    
    container = np.zeros((1, pop.shape[0]), str)
    indexes = [pop.shape[0]]
    expand_inputs = {}
    
    preprocess = {}
    level1 = {}
    level2 = {}
    level3 = {}
    correction = {}
    dic = 'preprocess'
    valid = ['preprocess_old', 'level1_old', 'level2_old', 'level3_old', 'correction_old']
    master = {}
    previous = ''
    counter = 0
    for key in map_genes:
        #key -= counter
        gene = map_genes[key]
        keys = list(gene.keys())
        node_name = re.search('([A-Za-z0-9]+)_', keys[0]).group(1)
        if not previous or 'end' in previous:
            previous = node_name
            
        if previous != node_name:
            for link in links[dic]:
                if previous in link[:len(previous)]:
                    connect = links[dic][link]
                    if type(connect) == list:
                        group = re.search('([A-Za-z]+)_', connect[0]).group(1)
                        vals = vars()[dic][group][connect[0]]
                        check = vars()[dic][group][connect[1]]
                        #WILL LIKELY NEED TO ALTER THIS AS WELL
                        #
                        #
                        #
                        #
                        vars()[dic][previous][link] = [val if check[c] else 0 for c, val in enumerate(vals)]
                    else:
                        try:
                            group = re.search('([A-Za-z]+)_', connect).group(1)
                            if group in vars()[dic]:
                                vars()[dic][previous][link] = vars()[dic][group][connect]
                            else:
                                for opt in valid:
                                    if group in vars()[opt]:
                                        break
                                vars()[dic][previous][link] = vars()[opt][group][connect]
                        except:
                            A=3
            if 'F' == previous[0]:
                expand_inputs[previous[1:]] = vars()[dic][previous]
                
            previous = node_name
            
        if 'end' in keys[0]:
            #NOTE: should be able to replace pholder with vars()[dic] and remove below two lines -> not working in spyder debug
            container, pholder, indexes = define_paths(container, vars()[dic], indexes)
            vars()[dic+'_old'] = vars()[dic].copy()
            vars()[dic].clear()
            vars()[dic].update(pholder)
            master[dic] = vars()[dic]
            dic = gene[keys[0]]#re.search('[A-Za-z]+_([A-Za-z]+)', keys[0]).group(1)
            counter += 1
            continue
        
        if node_name not in vars()[dic]:
            vars()[dic][node_name] = {}
            
        values = params[key-counter,:]
        
        isint = False
        if values.dtype == float:
            check_vals = [val.is_integer() for val in values]
            if sum(check_vals) == len(check_vals):
                isint = True
        #IF REST IN TASK SKIP CERTAIN VALUES
        const_ = np.all(values[0] == values)
        
        if keys[0][0] == '!':
            continue
        
        if len(gene) > 1 or '~construct~' in keys[0]:
            for l, i in enumerate(values):
                if round(i) in gene:
                    mapped = gene[round(i)]
                else:
                    mapped = i
                if type(mapped) == dict and keys[0][-1] == '_': #corresponds to mutually exclusive parameters
                    if False:#const_: #SHOULD BE ABLE TO DELETE THIS, IF USING UNDEFINED ABOVE
                        const_keys = list(mapped.keys())
                        const_values = mapped.values()
                        for i_, item_ in enumerate(const_values):
                            if item_: #depends on the fact that mutually exclusive will correspond to some non-zero or null value
                                mapped = {const_keys[i_]: item_}
                                break
                    if const_:
                        A=3
                    for k in mapped:
                        if keys[0][-1] == '_':
                            param = keys[0] + k
                        else:
                            param = node_name + '_' + k
                        if param not in vars()[dic][node_name]:
                             vars()[dic][node_name][param] = []
                             
                        vars()[dic][node_name][param].append(mapped[k])
                elif '~construct~' in keys[0]:
                    #ENUMERATE for access key to correct dictionary, add key
                    var_name = re.search('_([A-Za-z]+)', keys[0]).group(1)
                    key_name = re.search('_([A-Za-z]+)$', keys[0]).group(1)
                    param = node_name + '_' + var_name
                    
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    
                    if not isint:
                        rand = random.Random(i)
                        if type(mapped) == list:
                            if type(mapped[0]) == tuple or type(mapped[0]) == list:
                                mapped = [(m[0], rand.randint(m[1]-edge1, m[1]+edge2)) if len(m) == 2 else (m[0], rand.randint(m[1], m[2])) for m in mapped]
                            else:
                                print('UNDEFINED USE CASE')
                        elif mapped == tuple:
                            mapped = rand.randint(mapped[0], mapped[1])
                        else:
                            mapped = i
                    
                    if len(vars()[dic][node_name][param]) != pop.shape[0]:
                        construction_key = key_name
                        vars()[dic][node_name][param].append({key_name: mapped})
                    else:
                        if key_name in links:
                            if vars()[dic][node_name][param][l][construction_key] not in links[key_name]:
                                continue
                            
                        vars()[dic][node_name][param][l][key_name] = mapped
                else:
                    param = keys[0]
                        
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    
                    vars()[dic][node_name][param].append(mapped)
        else:
            if values.dtype == float:
                if isint:
                    values = values.astype(int)
                    
            vars()[dic][node_name][keys[0]] = values
    
    pipeline = generation * pop.shape[0]
    layout = BIDSLayout(data_dir)
    tasks = layout.get_tasks()
    for task in tasks:
        subjects = layout.get_subjects(task=task)[0]#:10]
        types = layout.get_datatypes()
        sessions = layout.get_sessions(task=task)
        runs = layout.get_runs(task=task)

        if 'anat' in types and 'func' in types:
            pipelines = analysis(exp_dir, working_dir, data_dir, out_dir)
            pipelines = pipelines.construct(subjects, sessions, runs, task, pipeline, master, expand_inputs)
            pipelines.inputs.inputnode.mask = mask
            pipelines.inputs.inputnode.task = task
            pipelines.inputs.inputnode.max_groups = group_num
            pipelines.run()
    
    A = 3
                            
                    

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
               
               gene_type=[float, 2],
               
               num_genes=num_genes,
               gene_space=gene_space,
               parent_selection_type='sss',
               keep_parents=-1,
               crossover_type="single_point",
               mutation_type="random",#look into adaptive
               mutation_probability=0.2,
               )
    ga.run()
    
if __name__ == "__main__":
    main()

