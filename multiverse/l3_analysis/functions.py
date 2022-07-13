#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:06:18 2021

@author: grahamseasons
"""
from functions import insert
from nipype.utils.functions import getsource
from workflows import write_out
import re

def group_contrast(copes, varcopes):
    outcopes = []
    outvarcopes = []
    multiple = type(copes[0])
    if multiple == list:
        numcon = len(copes[0])
    else:
        return [copes], [varcopes]
        
    for i in range(numcon):
        subcopes = [sub[i] for sub in copes]
        subvarcopes = [sub[i] for sub in varcopes]
        outcopes.append(subcopes)
        outvarcopes.append(subvarcopes)
        
    return outcopes, outvarcopes

def remove_container(cont):
    if isinstance(cont, list) and isinstance(cont[0], list) and len(cont) == 1:
        return cont[0]
    else:
        return cont
    
def mniMask(mask):
    import os
    old = os.path.join(os.getenv('FSLDIR'), 'data/standard/MNI152_T1_2mm.nii.gz')
    #old = os.path.join(os.getenv('FSLDIR'), 'data/standard/MNI152_T1_2mm_brain.nii.gz')
    if mask == old:
        mask = os.path.join(os.getenv('FSLDIR'), 'data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz')
    
    return mask

def t_test(covariate, subjects, demean):
    import pandas as pd
    import numpy as np
    covariate = pd.read_table(covariate)
    covariates = covariate.set_index('participant_id')
    covariates = covariates.loc[['sub-' + sub for sub in subjects]]
    categories = covariates.columns
    groupcat = 'group'
    EVs = {}
    contrasts = []
    
    #IF NO GROUPS OR ANYTHING
    if len(categories) > 0:
        if groupcat not in categories:
            groupcat = categories[0]
            
        #ASSUMES unpaired t-test
        #THIS IS EXPECTING STRING CATEGORIES FOR GROUPS -> could probably eliminate need with inclusion of json file
        group = covariates.groupby(groupcat)
        num_groups = len(group.count())
        group_ids = (group.ngroup() + 1).to_list()
        encoded = group.ngroup().to_list()
        labels = covariates[groupcat].unique()
        contrast = []
        
        for i in range(num_groups):
            ev = [1 if val == i else 0 for val in encoded]
            EVs[labels[i]] = ev
            
            solo = (labels[i] + ' mean', 'T', [labels[i]], [1]) #-> I THINK THIS MIGHT BE AN F CONTRAST
            contrast = [(labels[i] + '-' + lab, 'T', [labels[i], lab], [1,-1]) if lab != labels[i] else solo for lab in labels]
            contrasts += contrast #contrasts.append(contrast)
            
        #NOTE: FOR THE NON-GROUP COVARIATES THEY ARE ADDED AS IS RIGHT NOW -> NO DEMEANING/ORTHOGONALIZATION
        cov = covariates.drop(groupcat, axis=1)
        cat = categories.drop(groupcat)
        
        #orthog = np.array(list(EVs.values())).T

        EV_safe = EVs.copy()
        
        if demean:
            for c in cat:
                labels = cov[c].unique()
                if len(EV_safe) >= 1:#MIGHT NEED TO RETURN THIS TO > 1
                    for key in EV_safe:
                        if type(labels[0]) == str:
                            group = cov.groupby(c)
                            encoded = group.ngroup().to_list()
                            encoded_mean = np.mean(encoded)
                            #orthog, vec = find_orth(orthog, np.array(encoded, ndmin=2).T)
                            #EVs[c] = vec.tolist()
                            EVs[c + '_' + key] = [encoded[i] - encoded_mean if val else 0 for i, val in enumerate(EV_safe[key])]
                        else:
                            reg = cov[c].to_list()
                            reg_mean = np.mean(reg)
                            EVs[c + '_' + key] = [reg[i] - reg_mean if val else 0 for i, val in enumerate(EV_safe[key])]
                            #orthog, vec = find_orth(orthog, np.array(cov[c].to_list(), ndmin=2).T)
                            #EVs[c] = vec.tolist()
    else:
        single_group = [1] * len(covariates.index)
        label = 'group_mean'
        EVs[label] = single_group
        group_ids = single_group
        contrasts.append([label, 'T', [label], [1]])
        
            
    return EVs, contrasts, group_ids


def get_sink(inputs, files):
    if type(files) != list:
        files = [files]
        
    func_str = getsource(write_out)
    ind = func_str.find('):')
    params = ', ' + ', '.join(inputs)
    func_str = insert(func_str, ind, params)
    
    search = re.search('\n(\n)(\s+)(setattr)', func_str)
    ind = search.start(1)
    
    block = '\n' + search.group(2) + (search.group(2) + search.group(2)).join(["if isinstance(vars()[out], str) and os.path.isdir(vars()[out]):\n",
        "for file in {files}:\n",
        "    file_out = vars()[out] + '/' + file\n",
        "    if os.path.isdir(file_out) or os.path.isfile(file_out):\n",
        "        setattr(sink.inputs, 'pipelines/' + task + '.@' + str(i) + file, file_out)\n"])
        
    func_str = insert(func_str, search.start(3), "\n    "+search.group(2))
    func_str = insert(func_str, search.start(3), "else:")
    func_str = insert(func_str, ind, block)
        
    return func_str.format(files=str(files))