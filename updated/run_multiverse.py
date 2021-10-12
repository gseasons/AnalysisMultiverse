#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:41:49 2021

@author: grahamseasons
"""
import pygad as pg
import numpy as np
import re, os
from os.path import join as opj
from updated.functions import define_paths
from updated.analysis_pipeline import analysis
from bids.layout import BIDSLayout

exp_dir = '/Volumes/NewVolume/sup_preprocess_test_join'
working_dir = 'working_dir'
data_dir = '/Volumes/NewVolume/super_agers'
out_dir = exp_dir + '/processed'
mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
group_num = 10

links = {'decision_mc_mean': 'mcflirt_mean_vol',
         'Fregistration_wmthresh': ['Fregistration_wmthresh', 'Fregistration_bbr'],
         'Fregistration_regbbr_interp': ['Fregistration_regpre_interp', 'Fregistration_bbr'],
         'Fregistration_regbbr_no_resample': ['Fregistration_regpre_no_resample', 'Fregistration_bbr'],
         'Fregistration_applywarp_interp': 'Fregistration_regpre_interp',
         'Fregistration_applywarp_no_resample': 'Fregistration_regpre_no_resample',
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
         #RUN SEQUENTIALLY
         {'bet_frac': [0.5]}, #{'low': 0.15, 'high': 0.7, 'step': 0.05}},
         {'bet_vertical_gradient': [0.2]},#{'low': -0.5, 'high': 0.5, 'step': 0.05}},
         {'bet_': [1],#[0, 1, 2],
              0: {'robust': True, 'reduce_bias': False, 'remove_eyes': False},
              1: {'robust': False, 'reduce_bias': True, 'remove_eyes': False},
              2: {'robust': False, 'reduce_bias': False, 'remove_eyes': True}},
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
         {'Fregistration_bbr': [1],#**
              0: False,
              1: True},
         {'Fregistration_warplater': [0, 1],#**
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
         {'Fregistration_wmthresh': {'low': 0.2, 'high': 0.8, 'step': 0.1}},#**
         {'Fsmooth_susan': [1],#**
              0: False,
              1: True},
         {'Fsmooth_fwhm': range(2, 13, 2)},#**
         #ART
         {'end_preprocess': 'level1'}
         ]

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
    
    previous = ''
    counter = 0
    for key in map_genes:
        key -= counter
        gene = map_genes[key]
        keys = list(gene.keys())
        node_name = re.search('([A-Za-z]+)', keys[0]).group(1)
        if not previous or 'end' in previous:
            previous = node_name
            
        if previous != node_name:
            for link in links:
                if previous in link[:len(previous)]:
                    connect = links[link]
                    if type(connect) == list:
                        group = re.search('([A-Za-z]+)_', connect[0]).group(1)
                        vals = vars()[dic][group][connect[0]]
                        check = vars()[dic][group][connect[1]]
                        vars()[dic][previous][link] = [val if check[c] else 0 for c, val in enumerate(vals)]
                    else:
                        vars()[dic][previous][link] =  vars()[dic][re.search('([A-Za-z]+)_', connect).group(1)][connect]
            if 'F' == previous[0]:
                expand_inputs[previous[1:]] = vars()[dic][previous]
                
            previous = node_name
            
        if 'end' in keys[0]:
            #NOTE: should be able to replace pholder with vars()[dic] and remove below two lines -> not working in spyder debug
            container, pholder, indexes = define_paths(container, vars()[dic], indexes)
            vars()[dic].clear()
            vars()[dic].update(pholder)
            dic = re.search('[A-Za-z]+_([A-Za-z])', keys[0]).group(1)
            counter += 1
            continue
        
        if node_name not in vars()[dic]:
            vars()[dic][node_name] = {}
            
        values = params[key,:]
        const_ = np.all(values[0] == values)
        if len(gene) > 1:
            for i in values:
                mapped = gene[i]
                if type(mapped) == dict: #corresponds to mutually exclusive parameters
                    if const_:
                        const_keys = list(mapped.keys())
                        const_values = mapped.values()
                        for i_, item_ in enumerate(const_values):
                            if item_: #depends on the fact that mutually exclusive will correspond to some non-zero or null value
                                mapped = {const_keys[i_]: item_}
                                break
                    for k in mapped:
                        param = node_name + '_' + k
                        if param not in vars()[dic][node_name]:
                             vars()[dic][node_name][param] = []
                             
                        vars()[dic][node_name][param].append(mapped[k])
                else:
                    param = keys[0]
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    
                    vars()[dic][node_name][param].append(mapped)
        else:
            if values.dtype == float:
                check_vals = [val.is_integer() for val in values]
                if sum(check_vals) == len(check_vals):
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
            pipelines = pipelines.construct(subjects, sessions, runs, task, pipeline, preprocess, expand_inputs)
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

