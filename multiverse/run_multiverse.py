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

import shutil

from nipype.utils.profiler import log_nodes_cb
#META GETTING DELETED TOO EARLY? HOPEFULLY JUST A FUNCTION OF ACCIDENTALLY NOT USING IPYTHON
#TO DO: COMMENTS (especially GUI), CONTRIBUTE PLUGIN_BASE FILE (ALLOWS FOR DELETING USED NODES, IS SAFE FOR IDENTITY)
#       Update so that we request whole nodes instead of CPUs on any number of nodes -> should limit failures

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
        is_params = load('reproducibility', 'generation_'+str(generation)+'.pkl')
        if type(is_params) != str:
            params_ = is_params
            pop_ = params_.transpose()
            pipeline = generation * pop_.shape[0]
            solution_start = pipeline
    else:
        save('reproducibility', 'generation_'+str(generation)+'.pkl', params_)
    
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
            
        if config['rerun'] or checkpoints:
            frame = ''
        else:
            frame = load('', task+'.pkl')
        #TODO: FIGURE OUT BATCH LATER -> ENSURE HAVE ENOUGH SPACE, GB_PER_PIPE INACCURATE FOR HOW THINGS WOULD RUN (only valid if can parallelize literally everything)
        gb_per_pipe = len(subjects) * (len(sessions) + 1) * (len(runs) + 1) * 0.83
        
        if (config['storage'] / gb_per_pipe) < pop_.shape[0]:
            batch_size = int(config['storage'] / gb_per_pipe)
            iterations = math.ceil(pop_.shape[0] / batch_size)
        else:
            batch_size = pop_.shape[0]
            iterations = 1
        
        last_batch = 0
        batch_size = 8
        iterations = 2
        
        for batch in range(iterations):
            if checkpoints:
                workflows = glob.glob('/scratch/processed/reproducibility/' + task + '_workflow*')
                last_batch = len(workflows) - 1
                if batch < last_batch:
                    continue
                
            if (batch+1) * batch_size < pop_.shape[0]:
                params = params_[:, batch*batch_size:(batch+1)*batch_size]
                pop = pop_[batch*batch_size:(batch+1)*batch_size,:]
            else:
                params = params_[:, batch*batch_size:]
                pop = pop_[batch*batch_size:,:]
            print(params.shape)
            print(pop.shape)
                
            pipeline = pipeline_ + batch * batch_size
            
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
                
                plugin_args = {}
                import ipyparallel as ipp
                A = ipp.Cluster(n=6)
                plugin_args = {'cluster_id': A.cluster_id, 'profile': A.profile}
                #A.wait_for_engines
                A = A.start_cluster_sync()
                #print(A.ids)
                config['processing'] = 'IPython'
                #exit()
                if config['processing'] == 'SLURM':
                    config['processing'] = 'IPython'
                    plugin_args = {'profile': profile}
                    #CHANGE SO SAVES TO REPRODUCIBILITY
                if checkpoints and batch == last_batch:
                    pipelines = load('reproducibility', task + '_workflow_' + str(batch) + '.pkl')
                else:
                    save('reproducibility', task + '_workflow_' + str(batch) + '.pkl', pipelines)
                    
                pipelines.run(plugin=config['processing'], plugin_args=plugin_args)#plugin='Linear', plugin_args={'status_callback': log_nodes_cb})#plugin=config, plugin_args=plugin_args)
                if not config['debug']:
                    shutil.rmtree(working_dir)
                
    if 'num_generations' not in config:
        if checkpoints:
            os.rename('/scratch/processed/reproducibility/checkpoints', '/scratch/processed/reproducibility/checkpoints_finished')
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
        sol_per_pop = 16#config['pipelines']
    
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

