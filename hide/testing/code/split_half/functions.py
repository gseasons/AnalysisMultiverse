#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:55:54 2021

@author: grahamseasons
"""
from nipype.utils.functions import getsource
from functions import insert
from workflows import write_out

def calc(R, P, out_frame):
    import numpy as np
    import os, re
    import pandas
    import pickle
    perfect = np.array([1, 1])
    pipe = np.array([R, P])
    score = np.linalg.norm(perfect - pipe)
    pipeline_ = re.search('_i_([0-9]+)', os.getcwd()).group(1)
    with open(out_frame, 'rb') as f:
        frame = pickle.load(f)
        
    frame['R'][pipeline_] = R
    frame['P'][pipeline_] = P
    frame['Score'][pipeline_] = 1 - score
    
    with open(out_frame, 'wb') as f:
        pickle.dump(frame, f)
    
    return np.savetxt(os.getcwd() + '/score.txt', 1 - score)

def compare(stats, mask):
    from nilearn.input_data import NiftiMasker
    from nilearn.image import math_img
    from nilearn.plotting import plot_img_comparison as plot
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    plt.interactive(False)
    masker = NiftiMasker(math_img('img > 0', img=mask))
    masker.fit()
    num_groups = int(len(stats)/2)
    out_stats = []
    
    for i in range(num_groups):
        out = plot(stats[2*i], stats[2*i+1], masker, plot_hist=False)
        out_stats += out
    
    out = np.mean(out_stats)
    
    return out, np.savetxt(os.getcwd() + '/repro.txt', out) #, out_stats

def predict(covariate, bold, mask):
    import re
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.multioutput import MultiOutputClassifier
    import pandas as pd
    import numpy as np
    import nibabel as nib
    import os
    covariate = pd.read_table(covariate)
    covariates = covariate.set_index('participant_id')
    categories = covariates.columns
    data = covariates.dtypes
    
    reg = covariates.copy()
    
    for label in categories:
        if data[label] == 'O':
            reg = reg.drop(label, axis=1)
        else:
            covariates = covariates.drop(label, axis=1)
            
    regress = False
    classi = False
    
    if reg.shape[1] < len(categories):
        regress = True
        
    if covariates.shape[1] < len(categories):
        classi = True
            
    
    #covariates = covariate.columns
    num_cov = len(covariates.columns)
    #MAY NEED TO CHANGE SO THAT ADDITIONAL COLUMNS IN PARTICIPANTS ARE USED AS INPUTS AND ONLY DISEASE IS PREDICTED
    pred = []
    #TEST SPEED ON OLD FORMAT AS WELL
    mask_img = nib.load(mask)
    mask_ind = np.nonzero(mask_img.get_fdata())
    
    t_size = [nib.load(f).shape[-1] for f_c in bold for f in f_c]
    
    size = (int(np.ceil(np.ravel(bold).shape[0])/2), np.shape(mask_ind)[1]*max(t_size))
    dat0 = np.zeros(size).astype(np.int16)
    dat1 = np.zeros(size).astype(np.int16)
    subs = []
    
    for i, file_cont in enumerate(bold):
        ind = 0
        for file in file_cont:
            sub = re.search('_subject_([0-9S]+)', file).group(1)
            subs.append(sub)
            vars()['dat'+i][ind, :np.shape(mask_ind)[1]*t_size[ind]] = nib.load(file).get_fdata()[mask_ind].reshape(1, -1)
            ind += 1
            
    if num_cov:
        if classi:
            clf_c = RandomForestClassifier(random_state=2021, n_estimators=50)
        if regress:
            clf_r = RandomForestRegressor(random_state=2021, n_estimators=50)
        #change so that just assigns every column -> change frame into array
# =============================================================================
#         for group in groups:
#             X1 = group[0]
#             X2 = group[1]
#             if classi:
#                 y1 = []
#                 y2 = []
#             if regress:
#                 y1_r = []
#                 y2_r = []
#             
#             ind1 = []
#             ind2 = []
#             
#             index = []
#             for subgroup in group:
#                 for s in subgroup:
#                     index = 'sub-' + s
#                     group_ind = [pos for pos, val in enumerate(subs) if val == s]
#                     if s in X1:
#                         ind1 += group_ind
#                         if classi:
#                             y1.append(covariates.loc[index].to_list()*len(group_ind))
#                         if regress:
#                             y1_r.append(reg.loc[index].to_list()*len(group_ind))
#                     elif s in X2:
#                         ind2 += group_ind
#                         if classi:
#                             y2.append(covariates.loc[index].to_list()*len(group_ind))
#                         if regress:
#                             y2_r.append(reg.loc[index].to_list()*len(group_ind))
# =============================================================================
            
            if classi:
                Yc = covariates.to_numpy()
# =============================================================================
#                 Y1 = np.ravel(np.char.array(y1)).reshape(len(ind1), -1)
#                 Y2 = np.ravel(np.char.array(y2)).reshape(len(ind1), -1)
# =============================================================================
            if regress:
                Yr = reg.to_numpy()
# =============================================================================
#                 Y1_r = np.ravel(np.array(y1_r)).reshape(len(ind1), -1)
#                 Y2_r = np.ravel(np.array(y2_r)).reshape(len(ind1), -1)
# =============================================================================
                
# =============================================================================
#             dat1 = dat[tuple(ind1),...]
#             dat2 = dat[tuple(ind2),...]
# =============================================================================
            
            if classi:
                clf_c = clf_c.fit(dat0, Yc)
                prediction = clf_c.predict(dat1)
                pred.append((prediction == Yc).sum() / (Yc.shape[0]*Yc.shape[0]))
                
                clf_c = clf_c.fit(dat1, Yc)
                prediction = clf_c.predict(dat0)
                pred.append((prediction == Yc).sum() / (Yc.shape[0]*Yc.shape[0]))
                
            if regress:
                clf_r = clf_r.fit(dat0, Yr)
                prediction = clf_r.predict(dat1)
                pred.append((prediction == Yr).sum() / (Yr.shape[0]*Yr.shape[0]))
                
                clf_r = clf_r.fit(dat1, Yr)
                prediction = clf_r.predict(dat0)
                pred.append((prediction == Yr).sum() / (Yr.shape[0]*Yr.shape[0]))


            #clf = clf.fit(dat2, Y2)
            #pred.append(clf.score(dat1, Y1))
    else:
        print("NO PARTICIPANTS FILE")
        
    out = np.mean(pred)
        
    return out, np.savetxt(os.getcwd() + '/pred.txt', out) #, pred

def get_sink(inputs):
    func_str = getsource(write_out)
    ind = func_str.find('):')
    params = ', ' + ', '.join(inputs)
    func_str = insert(func_str, ind, params)
    return func_str

