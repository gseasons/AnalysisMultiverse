#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:41:49 2021

@author: grahamseasons
"""
import pygad as pg
import numpy as np
import os
from os.path import join as opj
from analysis_pipeline import analysis
from bids.layout import BIDSLayout
from functions import generate_dictionaries
import pickle
from nipype import config as conf
import json
import sys
from pathlib import Path

exp_dir = '/scratch'
working_dir = 'working_dir'
data_dir = '/data'
out_dir = exp_dir + '/processed'
mask = opj(os.getenv('FSLDIR'), 'data/standard/MNI152_T1_2mm_brain.nii.gz')
dir = os.path.dirname(os.path.abspath(__file__))

def fix_links(dic):
    links = {}
    for key in dic:
        links[key] = {}
        for link in dic[key]:
            if 'verify' in link:
                links[key].update({link['verify']: link['values']})
            elif 'node_to_add' in link:
                if 'on_off' in link:
                    links[key].update({link['node_to_add']: [link['node_to_copy'], link['on_off']]})
                else: 
                    links[key].update({link['node_to_add']: link['node_to_copy']})
            elif 'node_to_edit' in link:
                links[key].update({link['node_to_edit']: [link['node_to_edit'], link['on_off'], link['switch']]})
    return links

with open(opj(dir, 'configuration', 'multiverse_configuration.pkl'), 'rb') as f:
    genes = pickle.load(f)
    
    break_ = False
    for gene in genes:
        if 'l1d_bases' in gene.keys():
            for key in gene:
                if isinstance(gene[key], dict) and 'custom' in gene[key].keys():
                    for new_key in gene[key]['custom']:
                        if 'FSLDIR' in gene[key]['custom'][new_key]:
                            gene[key]['custom'][new_key] = gene[key]['custom'][new_key].format(FSLDIR=os.getenv('FSLDIR'))
                            break_ = True
                            break
                if break_:
                    break
        if break_:
            break
    
with open(opj(dir, 'configuration', 'general_configuration.pkl'), 'rb') as f:
    config = pickle.load(f)
    
with open(opj(dir, 'configuration', 'default_links.json')) as f:
    prelim_links = json.load(f)
    
links = fix_links(prelim_links)

conf.set("execution", "hash_method", "content")
conf.set("execution", "remove_node_directories", "false")
if not config['debug']:
    conf.set("execution", "remove_node_directories", "false")
    
    
wiggle = 10
map_genes = {}
solution_start = 0

layout = BIDSLayout(data_dir)
tasks = layout.get_tasks()

def load(path, file):
    out = os.path.join(out_dir, path, file)
    if os.path.isfile(out):
        with open(out, 'rb') as f:
            loaded = pickle.load(f)
    else:
        loaded = ''
    
    return loaded

def save(path, file, frame):
    out = os.path.join(out_dir, path)
    Path(out).mkdir(parents=True, exist_ok=True)
    out = os.path.join(out, file)
# =============================================================================
#     if not os.path.isdir(out_dir):
#         os.mkdir(out_dir)
# =============================================================================
        
    with open(out, 'wb') as f:
        pickle.dump(frame, f)
    
    return out

def check_pipes():
    unique = []
    for task in tasks:
        try:
            frame = load('', task+'.pkl')
            unique.append(frame.astype(str).drop_duplicates(subset=frame.columns.difference(['R', 'P', 'Score'])).shape[0] > config['pipelines'])
        except AttributeError:
            return
        
    if sum(unique) == config['pipelines']:
        return "stop"
    else:
        return
            

def fitness_func(solution, solution_idx):
    if config['split_half']:
        avg = []
        for task in tasks:
            frame = load('', task+'.pkl')
            avg.append(frame['Score'][solution_start+solution_idx])
        return np.mean(avg)
    else:
        return 1
    
def on_pop_gen(ga):
    gen = ga.generations_completed
    generation = gen
    pop = ga.population
    params = pop.transpose()
    pipeline = generation * pop.shape[0]
    global solution_start
    solution_start = pipeline
    
    if check_pipes():
        return "stop"
    
    if True:#config['rerun']:
        is_params = load('reproducibility', 'generation_'+str(generation)+'.pkl')
        if type(is_params) != str:
            params = is_params
            pop = params.transpose()
            pipeline = generation * pop.shape[0]
            solution_start = pipeline
    else:
        save('reproducibility', 'generation_'+str(generation)+'.pkl', params)
    
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
            
        if True:#config['rerun']:
            frame = ''
        else:
            frame = load('', task+'.pkl')
        
        master, expand_inputs, unique_pipelines = generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle, pipeline, frame)

        un = unique_pipelines.shape[0]
        test_unique = unique_pipelines.astype(str).drop_duplicates(subset=unique_pipelines.columns.difference(['R', 'P', 'Score']))
        test_un = test_unique.shape[0]
        
        if test_un < un:
            duplicates = unique_pipelines[unique_pipelines.astype(str).duplicated(keep=False, subset=unique_pipelines.columns.difference(['R', 'P', 'Score']))].astype(str)
            duplicates = duplicates.groupby(list(duplicates)).apply(lambda x: tuple(x.index)).to_list()
            for dup in duplicates:
                for row in dup:
                    if row == dup[0]:
                        continue
                    else:
                        unique_pipelines['R'][row] = unique_pipelines['R'][dup[0]]
                        unique_pipelines['P'][row] = unique_pipelines['P'][dup[0]]
                        unique_pipelines['Score'][row] = unique_pipelines['Score'][dup[0]]
                        
        to_run = [i for i in list(test_unique.index.values) if i >= pipeline]
        out_frame = save('', task+'.pkl', unique_pipelines)
        
        to_replace = [l-pipeline for l in range(pipeline, pipeline+pop.shape[0]) if l not in to_run]
        if to_replace:
            start_ind = min(to_run)
            params[:,to_replace] = params[:,start_ind-pipeline].reshape(-1,1)
            if start_ind != pipeline:
                params = np.delete(params, range(start_ind-pipeline), 1)
            pop = params.transpose()
            master, expand_inputs, _ = generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle, start_ind, frame)
        
        if 'anat' in types and 'func' in types and to_run:
            pipelines = analysis(exp_dir, working_dir, data_dir, out_dir)
            pipelines = pipelines.construct(subjects, sessions, runs, task, pipeline, master, expand_inputs, config['split_half'], to_run, config['networks'], out_frame)
            pipelines.inputs.inputnode.mask = mask
            pipelines.inputs.inputnode.task = task
            if config['processing'] == 'SLURM':
                config['processing'] = 'MultiProc'
            pipelines.run(plugin=config['processing'])
            
    if 'num_generations' not in config:
        sys.exit()
        #RUN PIPELINES IN BATCHES BASED ON NUMBER OF SUBJECTS/NUMBER OF PIPELINES -> max of ~0.83 GB PER SUBJECT PER PIPELINE
        #GA NOT SELECTED -> 1 generation, number of pipelines -> run in batches
        #ALLOW USER TO SPECIFY MAXIMUM OUTPUT FOLDER SIZE
        #TEST!        
        
        #TRY RUNNING ON COMPUTE CANADA

def main():
    if 'num_generations' in config:
        num_generations = config['num_generations']
        num_parents_mating = config['num_parents_mating']
        parent_selection_type = config['parent_selection_type']
        crossover_type = config['crossover_type']
        mutation_type = config['mutation_type']
        sol_per_pop = config['sol_per_pop']
    else:
        num_generations = 1
        num_parents_mating = 2
        parent_selection_type = 'random'
        crossover_type = 'single_point'
        mutation_type = 'random'
        sol_per_pop = config['pipelines']
    
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
               sol_per_pop=sol_per_pop,
               
               gene_type=[float, 3],
               
               num_genes=num_genes,
               gene_space=gene_space,
               parent_selection_type=parent_selection_type,
               
               keep_parents=0,
               save_solutions=False,
               
               crossover_type=crossover_type,
               mutation_type=mutation_type,
               mutation_probability=0.2,
               )
    ga.run()
    
if __name__ == "__main__":
    main()

