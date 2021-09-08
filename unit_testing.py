#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:48:57 2021

@author: grahamseasons
"""

def group_construction(subjects, split_half=True, average=False):
        import itertools, random
        group_container = []
        if split_half:
            prelim = list(itertools.combinations(subjects, round(len(subjects)/2)))
            pre_len = len(prelim)
            for i in range(pre_len):
                if not (pre_len % 2):
                    if i == (pre_len / 2):
                        break
                    else:
                        group_container.append([list(prelim[i]), list(prelim[-(i+1)])])
                else:
                    missed = [sub for sub in subjects if sub not in list(prelim[i])]
                    group_container.append([missed, list(prelim[i])])
        elif average:
            if type(subjects) == list:
                group_container.append([subjects])
            else:
                group_container.append([[subjects]])
        else:
            #PLACEHOLDER
            group_container.append([subjects, ['Single Group']])
            
        return random.sample(group_container, round(len(group_container)*0.1))


def construction(groups, copes, varcopes):
            import re
            subj_proc = len(copes)
            contrasts = len(copes[0])
            merge_contain_c = []
            merge_contain_v = []
            num_copes = []
            
            for group in groups:
                if group[-1][0] == 'Single Group':
                    for i in range(contrasts):
                        merge_contain_c.append([cope[i] for cope in copes])
                        merge_contain_v.append([varcope[i] for varcope in varcopes])
                        num_copes.append(len(group[0]))
                else:
                    for subset in group:
                        cont_c = []
                        cont_v = []
                        for i in range(contrasts):
                            cont_c.append([cope[i] for cope in copes if re.search('_subject_([0-9]+)', cope[i]).group(1) in subset])
                            cont_v.append([varcope[i] for varcope in varcopes if re.search('_subject_([0-9]+)', varcope[i]).group(1) in subset])
                            num_copes.append(len(subset))
                        
                        merge_contain_c.append(cont_c)
                        merge_contain_v.append(cont_v)       
            
            return merge_contain_c, merge_contain_v, num_copes
        
def t_test(covariate, subjects):
            import pandas as pd
            #PROBABLY PASS IN COVARIATE AS FILE NAME
            covariate = pd.read_table(covariate)
            covariates = covariate.set_index('participant_id')
            covariates = covariates.loc[['sub-' + sub for sub in subjects]]
            categories = covariates.columns
            groupcat = 'groups'
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
                    
                    solo = [labels[i] + ' mean', 'T', [labels[i]], [1]] #-> I THINK THIS MIGHT BE AN F CONTRAST
                    contrast = [[labels[i] + '-' + lab, 'T', [labels[i], lab], [1,-1]] if lab != labels[i] else solo for lab in labels]
                    contrasts.append(contrast)
                    
                #NOTE: FOR THE NON-GROUP COVARIATES THEY ARE ADDED AS IS RIGHT NOW -> NO DEMEANING/ORTHOGONALIZATION
                cov = covariates.drop(groupcat, axis=1)
                cat = categories.drop(groupcat)
                
                for c in cat:
                    labels = cov[c].unique()
                    
                    if type(labels[0]) == str:
                        encode = labels.ngroup().to_list()
                        EVs[c] = encode
                    else:
                        EVs[c] = cov[c].to_list()
            else:
                single_group = [1] * len(covariates.index)
                label = 'group_mean'
                EVs[label] = single_group
                group_ids = single_group
                contrasts.append([label, 'T', [label], [1]])
                
                    
            return EVs, contrasts, group_ids
        

def predict(groups, covariate, bold, mask):
            import re
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.multioutput import MultiOutputClassifier
            import pandas as pd
            import numpy as np
            import nibabel as nib
            import time
            start = time.time()
            covariate = pd.read_table(covariate)
            covariates = covariate.set_index('participant_id')
            #covariates = covariate.columns
            num_cov = len(covariates.columns)
            
            pred = []
            #TEST SPEED ON OLD FORMAT AS WELL
            mask_img = nib.load(mask)
            mask_ind = np.nonzero(mask_img.get_fdata())
            
            t_size = nib.load(bold[0][0]).shape[-1]
            
            size = (int(np.ceil(np.ravel(bold).shape[0])), np.shape(mask_ind)[1]*t_size)
            dat = np.zeros(size)
            subs = []
            ind = 0
            for file_cont in bold:
                for file in file_cont:
                    sub = re.search('_subject_([0-9]+)', file).group(1)
                    subs.append(sub)
                    dat[ind, :] = nib.load(file).get_fdata()[mask_ind].reshape(1, -1)
                    ind += 1
                    
            if num_cov:
                svm = RandomForestClassifier(random_state=2021, n_estimators=50)
                if num_cov > 1:
                    clf = MultiOutputClassifier(svm)
                else:
                    clf = svm
            #IMPROVE SPEED WITH INDEXING -> READ IN EVERYTHING ONLY ONCE
                for group in groups:
                    X1 = group[0]
                    X2 = group[1]
                    y1 = []
                    y2 = []
                    
                    ind1 = []
                    ind2 = []
                    
                    index = []
                    for subgroup in group:
                        for s in subgroup:
                            index = 'sub-' + s
                            group_ind = [pos for pos, val in enumerate(subs) if val == s]
                            if s in X1:
                                ind1 += group_ind
                                y1.append(covariates.loc[index].to_list()*len(group_ind))
                            elif s in X2:
                                ind2 += group_ind
                                y2.append(covariates.loc[index].to_list()*len(group_ind))
                            
                    Y1 = np.ravel(np.char.array(y1))
                    Y2 = np.ravel(np.char.array(y2))
                    
                    dat1 = dat[tuple(ind1),...]#np.take(dat, ind1, axis=0)
                    dat2 = dat[tuple(ind2),...]#np.take(dat, ind2, axis=0)
                    
                    clf = clf.fit(dat1, Y1)
                    pred.append(clf.score(dat2, Y2))
                    
                    clf = clf.fit(dat2, Y2)
                    pred.append(clf.score(dat1, Y1))
                    print(time.time() - start)
            else:
                print("NO PARTICIPANTS FILE")
                
            return pred, np.mean(pred)
        
def compare(stats, mask):
            from nilearn.input_data import NiftiMasker
            from nilearn.image import math_img
            from nilearn.plotting import plot_img_comparison as plot
            import numpy as np
            import matplotlib.pyplot as plt
            plt.interactive(False)
            masker = NiftiMasker(math_img('img > 0', img=mask))
            masker.fit()
            num_groups = int(len(stats)/2)
            out_stats = []
            
            for i in range(num_groups):
                out = plot(stats[2*i], stats[2*i+1], masker, plot_hist=False)
                out_stats += out
                break
            
            return out_stats, np.mean(out_stats)
        
from os.path import join as opj
import os

bold = [['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_04/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_04/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_05/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_05/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_06/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_06/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_07/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_07/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_08/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_08/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_09/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_09/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz']]
#[['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz']]
covariate = '/Users/grahamseasons/fMRI/test_data/ds000114/participants.tsv'
groups = [[['01', '02', '03', '04', '05'], ['06', '07', '08', '09', '10']], [['01', '02', '03', '04', '06'], ['05', '07', '08', '09', '10']], [['01', '02', '03', '04', '07'], ['05', '06', '08', '09', '10']], [['01', '02', '03', '04', '08'], ['05', '06', '07', '09', '10']], [['01', '02', '03', '04', '09'], ['05', '06', '07', '08', '10']], [['01', '02', '03', '04', '10'], ['05', '06', '07', '08', '09']], [['01', '02', '03', '05', '06'], ['04', '07', '08', '09', '10']], [['01', '02', '03', '05', '07'], ['04', '06', '08', '09', '10']], [['01', '02', '03', '05', '08'], ['04', '06', '07', '09', '10']], [['01', '02', '03', '05', '09'], ['04', '06', '07', '08', '10']], [['01', '02', '03', '05', '10'], ['04', '06', '07', '08', '09']], [['01', '02', '03', '06', '07'], ['04', '05', '08', '09', '10']], [['01', '02', '03', '06', '08'], ['04', '05', '07', '09', '10']], [['01', '02', '03', '06', '09'], ['04', '05', '07', '08', '10']], [['01', '02', '03', '06', '10'], ['04', '05', '07', '08', '09']], [['01', '02', '03', '07', '08'], ['04', '05', '06', '09', '10']], [['01', '02', '03', '07', '09'], ['04', '05', '06', '08', '10']], [['01', '02', '03', '07', '10'], ['04', '05', '06', '08', '09']], [['01', '02', '03', '08', '09'], ['04', '05', '06', '07', '10']], [['01', '02', '03', '08', '10'], ['04', '05', '06', '07', '09']], [['01', '02', '03', '09', '10'], ['04', '05', '06', '07', '08']], [['01', '02', '04', '05', '06'], ['03', '07', '08', '09', '10']], [['01', '02', '04', '05', '07'], ['03', '06', '08', '09', '10']], [['01', '02', '04', '05', '08'], ['03', '06', '07', '09', '10']], [['01', '02', '04', '05', '09'], ['03', '06', '07', '08', '10']], [['01', '02', '04', '05', '10'], ['03', '06', '07', '08', '09']], [['01', '02', '04', '06', '07'], ['03', '05', '08', '09', '10']], [['01', '02', '04', '06', '08'], ['03', '05', '07', '09', '10']], [['01', '02', '04', '06', '09'], ['03', '05', '07', '08', '10']], [['01', '02', '04', '06', '10'], ['03', '05', '07', '08', '09']], [['01', '02', '04', '07', '08'], ['03', '05', '06', '09', '10']], [['01', '02', '04', '07', '09'], ['03', '05', '06', '08', '10']], [['01', '02', '04', '07', '10'], ['03', '05', '06', '08', '09']], [['01', '02', '04', '08', '09'], ['03', '05', '06', '07', '10']], [['01', '02', '04', '08', '10'], ['03', '05', '06', '07', '09']], [['01', '02', '04', '09', '10'], ['03', '05', '06', '07', '08']], [['01', '02', '05', '06', '07'], ['03', '04', '08', '09', '10']], [['01', '02', '05', '06', '08'], ['03', '04', '07', '09', '10']], [['01', '02', '05', '06', '09'], ['03', '04', '07', '08', '10']], [['01', '02', '05', '06', '10'], ['03', '04', '07', '08', '09']], [['01', '02', '05', '07', '08'], ['03', '04', '06', '09', '10']], [['01', '02', '05', '07', '09'], ['03', '04', '06', '08', '10']], [['01', '02', '05', '07', '10'], ['03', '04', '06', '08', '09']], [['01', '02', '05', '08', '09'], ['03', '04', '06', '07', '10']], [['01', '02', '05', '08', '10'], ['03', '04', '06', '07', '09']], [['01', '02', '05', '09', '10'], ['03', '04', '06', '07', '08']], [['01', '02', '06', '07', '08'], ['03', '04', '05', '09', '10']], [['01', '02', '06', '07', '09'], ['03', '04', '05', '08', '10']], [['01', '02', '06', '07', '10'], ['03', '04', '05', '08', '09']], [['01', '02', '06', '08', '09'], ['03', '04', '05', '07', '10']], [['01', '02', '06', '08', '10'], ['03', '04', '05', '07', '09']], [['01', '02', '06', '09', '10'], ['03', '04', '05', '07', '08']], [['01', '02', '07', '08', '09'], ['03', '04', '05', '06', '10']], [['01', '02', '07', '08', '10'], ['03', '04', '05', '06', '09']], [['01', '02', '07', '09', '10'], ['03', '04', '05', '06', '08']], [['01', '02', '08', '09', '10'], ['03', '04', '05', '06', '07']], [['01', '03', '04', '05', '06'], ['02', '07', '08', '09', '10']], [['01', '03', '04', '05', '07'], ['02', '06', '08', '09', '10']], [['01', '03', '04', '05', '08'], ['02', '06', '07', '09', '10']], [['01', '03', '04', '05', '09'], ['02', '06', '07', '08', '10']], [['01', '03', '04', '05', '10'], ['02', '06', '07', '08', '09']], [['01', '03', '04', '06', '07'], ['02', '05', '08', '09', '10']], [['01', '03', '04', '06', '08'], ['02', '05', '07', '09', '10']], [['01', '03', '04', '06', '09'], ['02', '05', '07', '08', '10']], [['01', '03', '04', '06', '10'], ['02', '05', '07', '08', '09']], [['01', '03', '04', '07', '08'], ['02', '05', '06', '09', '10']], [['01', '03', '04', '07', '09'], ['02', '05', '06', '08', '10']], [['01', '03', '04', '07', '10'], ['02', '05', '06', '08', '09']], [['01', '03', '04', '08', '09'], ['02', '05', '06', '07', '10']], [['01', '03', '04', '08', '10'], ['02', '05', '06', '07', '09']], [['01', '03', '04', '09', '10'], ['02', '05', '06', '07', '08']], [['01', '03', '05', '06', '07'], ['02', '04', '08', '09', '10']], [['01', '03', '05', '06', '08'], ['02', '04', '07', '09', '10']], [['01', '03', '05', '06', '09'], ['02', '04', '07', '08', '10']], [['01', '03', '05', '06', '10'], ['02', '04', '07', '08', '09']], [['01', '03', '05', '07', '08'], ['02', '04', '06', '09', '10']], [['01', '03', '05', '07', '09'], ['02', '04', '06', '08', '10']], [['01', '03', '05', '07', '10'], ['02', '04', '06', '08', '09']], [['01', '03', '05', '08', '09'], ['02', '04', '06', '07', '10']], [['01', '03', '05', '08', '10'], ['02', '04', '06', '07', '09']], [['01', '03', '05', '09', '10'], ['02', '04', '06', '07', '08']], [['01', '03', '06', '07', '08'], ['02', '04', '05', '09', '10']], [['01', '03', '06', '07', '09'], ['02', '04', '05', '08', '10']], [['01', '03', '06', '07', '10'], ['02', '04', '05', '08', '09']], [['01', '03', '06', '08', '09'], ['02', '04', '05', '07', '10']], [['01', '03', '06', '08', '10'], ['02', '04', '05', '07', '09']], [['01', '03', '06', '09', '10'], ['02', '04', '05', '07', '08']], [['01', '03', '07', '08', '09'], ['02', '04', '05', '06', '10']], [['01', '03', '07', '08', '10'], ['02', '04', '05', '06', '09']], [['01', '03', '07', '09', '10'], ['02', '04', '05', '06', '08']], [['01', '03', '08', '09', '10'], ['02', '04', '05', '06', '07']], [['01', '04', '05', '06', '07'], ['02', '03', '08', '09', '10']], [['01', '04', '05', '06', '08'], ['02', '03', '07', '09', '10']], [['01', '04', '05', '06', '09'], ['02', '03', '07', '08', '10']], [['01', '04', '05', '06', '10'], ['02', '03', '07', '08', '09']], [['01', '04', '05', '07', '08'], ['02', '03', '06', '09', '10']], [['01', '04', '05', '07', '09'], ['02', '03', '06', '08', '10']], [['01', '04', '05', '07', '10'], ['02', '03', '06', '08', '09']], [['01', '04', '05', '08', '09'], ['02', '03', '06', '07', '10']], [['01', '04', '05', '08', '10'], ['02', '03', '06', '07', '09']], [['01', '04', '05', '09', '10'], ['02', '03', '06', '07', '08']], [['01', '04', '06', '07', '08'], ['02', '03', '05', '09', '10']], [['01', '04', '06', '07', '09'], ['02', '03', '05', '08', '10']], [['01', '04', '06', '07', '10'], ['02', '03', '05', '08', '09']], [['01', '04', '06', '08', '09'], ['02', '03', '05', '07', '10']], [['01', '04', '06', '08', '10'], ['02', '03', '05', '07', '09']], [['01', '04', '06', '09', '10'], ['02', '03', '05', '07', '08']], [['01', '04', '07', '08', '09'], ['02', '03', '05', '06', '10']], [['01', '04', '07', '08', '10'], ['02', '03', '05', '06', '09']], [['01', '04', '07', '09', '10'], ['02', '03', '05', '06', '08']], [['01', '04', '08', '09', '10'], ['02', '03', '05', '06', '07']], [['01', '05', '06', '07', '08'], ['02', '03', '04', '09', '10']], [['01', '05', '06', '07', '09'], ['02', '03', '04', '08', '10']], [['01', '05', '06', '07', '10'], ['02', '03', '04', '08', '09']], [['01', '05', '06', '08', '09'], ['02', '03', '04', '07', '10']], [['01', '05', '06', '08', '10'], ['02', '03', '04', '07', '09']], [['01', '05', '06', '09', '10'], ['02', '03', '04', '07', '08']], [['01', '05', '07', '08', '09'], ['02', '03', '04', '06', '10']], [['01', '05', '07', '08', '10'], ['02', '03', '04', '06', '09']], [['01', '05', '07', '09', '10'], ['02', '03', '04', '06', '08']], [['01', '05', '08', '09', '10'], ['02', '03', '04', '06', '07']], [['01', '06', '07', '08', '09'], ['02', '03', '04', '05', '10']], [['01', '06', '07', '08', '10'], ['02', '03', '04', '05', '09']], [['01', '06', '07', '09', '10'], ['02', '03', '04', '05', '08']], [['01', '06', '08', '09', '10'], ['02', '03', '04', '05', '07']], [['01', '07', '08', '09', '10'], ['02', '03', '04', '05', '06']]]#[[['01', '02'], ['03', '10']], [['01', '03'], ['02', '10']], [['01', '10'], ['02', '03']]]
mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
predict(groups, covariate, bold, mask)

stats = [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo0/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo1/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo2/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo3/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo4/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo5/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo6/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo7/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo8/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo9/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo10/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo11/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo12/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo13/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo14/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo15/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo16/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo17/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo18/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo19/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo20/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo21/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo22/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo23/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo24/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo25/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo26/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo27/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo28/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo29/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo30/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo31/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo32/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo33/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo34/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo35/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo36/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo37/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo38/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo39/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo40/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo41/stats/zstat1.nii.gz']]
#compare(stats, mask)    
#copes = [[['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']]]
#groups = [[['01', '02'], ['03', '10']], [['01', '03'], ['02', '10']], [['01', '10'], ['02', '03']]]
#varcopes = [[['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']]]
#t_test('/Users/grahamseasons/fMRI/test_data/ds000114/participants.tsv', ['01', '02', '03', '04'])
#construction(groups, copes, varcopes)
#group_construction(['01','02','03','04'])