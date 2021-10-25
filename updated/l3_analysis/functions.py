#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:06:18 2021

@author: grahamseasons
"""
def t_test(covariate, subjects):
    import pandas as pd
    #import numpy as np
    #from functions import find_orth
    #PROBABLY PASS IN COVARIATE AS FILE NAME
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
        #THIS IS EXPECTING STRING CATEGORIES FOR GROUPS -> could probably eliminate need with inclusing of json file
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
        
        for c in cat:
            labels = cov[c].unique()
            if len(EV_safe) > 1:
                for key in EV_safe:
                    if type(labels[0]) == str:
                        group = cov.groupby(c)
                        encoded = group.ngroup().to_list()
                        #orthog, vec = find_orth(orthog, np.array(encoded, ndmin=2).T)
                        #EVs[c] = vec.tolist()
                        EVs[c + '_' + key] = [encoded[i] if val else 0 for i, val in enumerate(EV_safe[key])]
                    else:
                        reg = cov[c].to_list()
                        EVs[c + '_' + key] = [reg[i] if val else 0 for i, val in enumerate(EV_safe[key])]
                        #orthog, vec = find_orth(orthog, np.array(cov[c].to_list(), ndmin=2).T)
                        #EVs[c] = vec.tolist()
    else:
        single_group = [1] * len(covariates.index)
        label = 'group_mean'
        EVs[label] = single_group
        group_ids = single_group
        contrasts.append([label, 'T', [label], [1]])
        
            
    return EVs, contrasts, group_ids