#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:41:49 2021

@author: grahamseasons
"""
import pygad as pg
import numpy as np
import os, math
from os.path import join as opj
from analysis_pipeline import analysis
from bids.layout import BIDSLayout
from functions import generate_dictionaries
import pickle
from nipype import config as conf
import json
import sys
import glob
from pathlib import Path
import pandas as pd
import re

import shutil

from nipype.utils.profiler import log_nodes_cb
#TO DO: COMMENTS (especially GUI), CONTRIBUTE PLUGIN_BASE FILE (ALLOWS FOR DELETING USED NODES, IS SAFE FOR IDENTITY)
#       Update so that we request whole nodes instead of CPUs on any number of nodes -> should limit failures
# ENABLE USER TO SELECT BATCH SIZES IN GUI (HOW MANY RUNS THEY WANT PIPELINES TO BE SPLIT ACROSS)

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
                            gene[key]['custom'][new_key] = "{FSLDIR}/etc/default_flobs.flobs".format(FSLDIR=os.getenv('FSLDIR'))#gene[key]['custom'][new_key].format(FSLDIR=os.getenv('FSLDIR'))
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
#conf.set("execution", "remove_node_directories", "false")#TRUE USUALLY


# =============================================================================
# conf.enable_resource_monitor()
# import logging
# callback_log_path = '/scratch/run_stats.log'
# logging.basicConfig(filename=callback_log_path, level=logging.DEBUG)
# logger = logging.getLogger('callback')
# handler = logging.FileHandler(callback_log_path)
# logger.addHandler(handler)
# =============================================================================

if not config['debug']:
    conf.set("execution", "remove_node_directories", "true")
    
if sys.argv[1] == "True":
    config['rerun'] = True
else:
    config['rerun'] = False
    
if len(sys.argv) > 2:
    profile = sys.argv[2]
    
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

def organize(task, out_frame):
    """Creates a dictionary of final output files, and parameters for each pipeline - excludes parameters that are unchanged across all pipelines
       
       Structure:
           {pipeline: {network: {contrast: file}},
                      {parameters: {parameters}}
                      }
    """
    processed = {'pipeline': {}}
    pathlist = Path(out_dir).glob('**/*_corrected_[0-9]*')
    dat_frame = out_frame
    
    with open(dat_frame, 'rb') as file:
        dat_frame = pickle.load(file)
        
    comp = pd.DataFrame(np.roll(dat_frame.values, 1, axis=0), index=dat_frame.index)
        
    for path in pathlist:
        path = str(path)
        network = int(re.search('.*_network_([0-9]+)', path).group(1))
        contrast = int(re.search('.*_corrected_([0-9]+).nii.gz', path).group(1))
        pipeline = int(re.search('.*_i_([0-9]+)', path).group(1))
        if pipeline in processed['pipeline']:
            if network in processed['pipeline'][pipeline]['network']:
                processed['pipeline'][pipeline]['network'][network]['contrast'][contrast] = path
            else:
                processed['pipeline'][pipeline]['network'][network] = {'contrast': {contrast: path}}
        else:
            processed['pipeline'][pipeline] = {'network': {network: {'contrast': {contrast: path}}}}
            
        pipe_dat = dat_frame.loc[pipeline]
        for i, column in enumerate(dat_frame):
            col = pipe_dat[column]
            if (comp[i] == dat_frame[column]).all():
                continue
            
            if 'parameters' not in processed['pipeline'][pipeline]:
                processed['pipeline'][pipeline]['parameters'] = {}
            
            if isinstance(col, dict):
                for key in col:
                    processed['pipeline'][pipeline]['parameters'][key] = col[key]
            else:
                processed['pipeline'][pipeline]['parameters'][column] = col
    
    return save('', task+'_organized.pkl', processed)

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
    pop_ = ga.population
    params_ = pop_.transpose()
    pipeline_ = generation * pop_.shape[0]
    global solution_start
    solution_start = pipeline_
    
    if check_pipes():
        return "stop"
    
    checkpoints = glob.glob('/scratch/processed/reproducibility/checkpoints/checkpoint_*.pkl')
    
    if config['rerun'] or checkpoints:
        generation = glob.glob('/scratch/processed/reproducibility/generation_*.pkl')
        generation = len(generation) - 1
        is_params = load('reproducibility', 'generation_'+str(generation)+'.pkl')
        if type(is_params) != str:
            params_ = is_params
            pop_ = params_.transpose()
            pipeline_ = generation * pop_.shape[0]
            solution_start = pipeline_
    else:
        save('reproducibility', 'generation_'+str(generation)+'.pkl', params_)
    
    for task in tasks:
        subjects = layout.get_subjects(task=task)
        subjects.sort()
        types = layout.get_datatypes()
        sessions = layout.get_sessions(task=task)
        runs = layout.get_runs(task=task)
        
        if sessions or runs:
            multiscan = True
        else:
            multiscan = False
            
        if config['rerun'] or checkpoints:
            frame = ''
        else:
            frame = load('', task+'.pkl')
        #TODO: FIGURE OUT BATCH LATER -> ENSURE HAVE ENOUGH SPACE, GB_PER_PIPE INACCURATE FOR HOW THINGS WOULD RUN (only valid if can parallelize literally everything)
        gb_per_pipe = len(subjects) * (len(sessions) + 1) * (len(runs) + 1) * 0.83
        
        batch_size = config['batches']
        iterations = math.ceil(pop_.shape[0] / batch_size)
        
        if (config['storage'] / gb_per_pipe) < pop_.shape[0]:
            batch_size_ = int(config['storage'] / gb_per_pipe)
            iterations_ = math.ceil(pop_.shape[0] / batch_size)
        else:
            batch_size_ = pop_.shape[0]
            iterations_ = 1
            
        if batch_size > batch_size_:
            batch_size = batch_size_
            iterations = iterations_
        
        for batch in range(iterations):
            if checkpoints:
                workflows = glob.glob('/scratch/processed/reproducibility/' + task + '_workflow*')
                last_batch = len(workflows) - 1
                if batch < last_batch:
                    continue
                elif batch == last_batch:
                    frame = load('', task+'.pkl')
                    frame.drop(range(batch*batch_size, (batch+1)*batch_size))
                
            if (batch+1) * batch_size < pop_.shape[0]:
                params = params_[:, batch*batch_size:(batch+1)*batch_size]
                pop = pop_[batch*batch_size:(batch+1)*batch_size,:]
            else:
                params = params_[:, batch*batch_size:]
                pop = pop_[batch*batch_size:,:]
                
            pipeline = pipeline_ + batch * batch_size
            
            if type(frame) != str:
                already_run = set(range(frame.shape[0]))
            else:
                already_run = set()
            master, expand_inputs, unique_pipelines = generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle, pipeline, frame)
            
            un = unique_pipelines.shape[0]
            test_unique = unique_pipelines.astype(str).drop_duplicates(subset=unique_pipelines.columns.difference(['R', 'P', 'Score']))
            test_un = test_unique.shape[0]
            
            if test_un < un:
                duplicates = unique_pipelines[unique_pipelines.astype(str).duplicated(keep=False, subset=unique_pipelines.columns.difference(['R', 'P', 'Score']))].astype(str)
                duplicates = duplicates.groupby(list(duplicates)).apply(lambda x: tuple(x.index)).to_list()
                for dup in duplicates:
                    for row in dup:
                        already_run.add(row)
                        
                        if row == dup[0]:
                            continue
                        else:
                            unique_pipelines['R'][row] = unique_pipelines['R'][dup[0]]
                            unique_pipelines['P'][row] = unique_pipelines['P'][dup[0]]
                            unique_pipelines['Score'][row] = unique_pipelines['Score'][dup[0]]
            
            to_run = [i for i in list(test_unique.index.values) if i >= pipeline and i not in already_run]
            
            if not to_run:
                continue
            
            frame = unique_pipelines
            
            to_replace = [l-pipeline for l in range(pipeline, pipeline+pop.shape[0]) if l not in to_run]
            
            if to_replace:
                start_ind = min(to_run)
                try:
                    params[:,to_replace] = params[:,start_ind-pipeline].reshape(-1,1)
                except IndexError:
                    continue
                if start_ind != pipeline:
                    params = np.delete(params, range(start_ind-pipeline), 1)
                pop = params.transpose()
                master, expand_inputs, _ = generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle, start_ind, frame)
            
            out_frame = save('', task+'.pkl', unique_pipelines)
            
            if 'anat' in types and 'func' in types and to_run:
                pipelines = analysis(exp_dir, working_dir, data_dir, out_dir)
                pipelines = pipelines.construct(subjects, sessions, runs, task, pipeline, master, expand_inputs, config['split_half'], to_run, config['networks'], out_frame)
                pipelines.inputs.inputnode.mask = mask
                pipelines.inputs.inputnode.task = task
                
                plugin_args = {}

                if config['processing'] == 'SLURM':
                    config['processing'] = 'IPython'
                    plugin_args = {'profile': profile}
                    
                if checkpoints and batch == last_batch:
                    pipelines = load('reproducibility', task + '_workflow_' + str(batch) + '.pkl')
                else:
                    save('reproducibility', task + '_workflow_' + str(batch) + '.pkl', pipelines)
                    
                pipelines.run(plugin=config['processing'], plugin_args=plugin_args)
                
                #USE IN FINAL DATA ANALYSIS
                organized = organize(task, out_frame)
                
                if not config['debug']:
                    if os.path.exists('/scratch/processed/reproducibility/checkpoints'):
                        os.rename('/scratch/processed/reproducibility/checkpoints', '/scratch/processed/reproducibility/checkpoints_' + task + '_batch_' + str(batch))
                    shutil.rmtree(working_dir)
                
    if 'num_generations' not in config:
        if checkpoints:
            os.rename('/scratch/processed/reproducibility/checkpoints', '/scratch/processed/reproducibility/checkpoints_finished')
        
        #TODO: DATA ANALYSIS
        sys.exit()#return "stop"
        #RUN PIPELINES IN BATCHES BASED ON NUMBER OF SUBJECTS/NUMBER OF PIPELINES -> max of ~0.83 GB PER SUBJECT PER PIPELINE
        #GA NOT SELECTED -> 1 generation, number of pipelines -> run in batches
        #ALLOW USER TO SPECIFY MAXIMUM OUTPUT FOLDER SIZE
        #TEST!        
        
        #TODO: TEST THAT DELETE FILES WHEN DONE WORKS WITH IPYTHON
        #      BENCHMARK WORKFLOW -> ~resources for each node
        #      REPLACE NIPYPE FILES WITH CUSTOM -> either in build or via bind? YES VIA BIND
        #      TRY RUNNING SMALL SET ON COMPUTE CANADA TOP TO BOTTOM
        #      MAKE SURE DOCKER IN BACKGROUND WORKS
        #      COMPILE SOURCES INTO ZOTERO
        #      DO RESOURCE PROFILING ON COMPUTE CANADA TO SEE WHY CPU EFFICIENCY SO LOW
        
    #*#*#FIX PATHS -> EACH IPYTHON ENGINE NEEDS TO BE LAUNCHED WITH THE SAME PATHS -> CURRENTLY GETTING SOME WHICH DON'T KNOW WHERE TO FIND USER DEFINED PACKAGES
        
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
    
    #CAN ADD POST PROCESSING THINGS HERE
    
if __name__ == "__main__":
    main()

