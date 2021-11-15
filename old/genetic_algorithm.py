#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 12:50:37 2021

@author: grahamseasons
"""
import pygad as pg
import numpy as np
frac_mask = [0.3, 0.3]
discard = [4, 4]
fwhm = [4, 8]
bet_frac = [0.5, 0.5]
robust = [True, True]
wm_thresh = [0.5, 0.5]
dof_f = [6, 6]
interp = ['spline', 'spline'] #0: 'spline', 1: 'nn', 2: 'sinc'
#iso = [4, 4]
bbr = [True, True]
susan = [True, True]
HP = [128, 140]
#serial_cor = [True, True]
#base_switch = [False, False]
#gamma = [False, False]
#base_val = [3, 3]#connected but unused
#mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
mode_l2 = ['fe', 'fe']
mode_l3 = ['flame1', 'flame1']
method = ['fwe', 'cluster']#'fdr']#'cluster']
p = [0.025, 0.025]
connectivity = [26, 26]
z_thresh = [3.2, 2.6]
warp_post_feat = [True, True]
generation = 0
preproc = {'FLIRT': {'free_param': []},
           'SMOOTHING': {'susan': [],#THESE WILL BE THE ACTUAL CONNECTION NAMES FOR EACH
                         'fwhm': []}
           }

genes = {'0': 'frac_mask',
         '1': 'discard',
         '2': 'fwhm',
         '3': 'bet_frac',
         '4': 'robust',
         '5': 'wm_thresh',
         '6': 'dof_f',
         '7': 'interp',
         '8': 'bbr',
         '9': 'susan',
         }


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
    preproc_act = preproc
    for i in range(params.shape[0]):
        if i == 2 or i == 9:#i < 10:
            preproc['SMOOTHING'][genes[str(i)]] = params[i,:]
    A=3
    #CONSTRUCT PIPELINES, 
    
    
    

def main():
    num_generations = 5
    num_parents_mating = 4
    fitness_function = fitness_func
    num_genes = 18
    gene_space = [{"low": 0.1, "high": 0.7},
                  [4, 5],
                  range(2, 13),
                  {"low": 0.2, "high": 0.7},
                  [0, 1],
                  {"low": 0.2, "high": 0.7},
                  [4, 6, 8],
                  [0, 1, 2],
                  [0, 1],
                  [0, 1],
                  range(10, 200),
                  [0, 1],
                  [0, 1, 2, 3],
                  [0, 1, 2, 3],
                  [0, 1, 2],
                  [18, 26, 32],
                  {"low": 2.3, "high": 3.2},
                  [0, 1]
                  ]
    allow_duplicate_genes = False
    ga = pg.GA(num_generations=num_generations,
                        num_parents_mating=num_parents_mating,
                        fitness_func=fitness_func,
                        on_start=on_pop_gen,
                        on_generation=on_pop_gen,
                        sol_per_pop=8,
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


















