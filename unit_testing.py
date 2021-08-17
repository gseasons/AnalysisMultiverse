#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:48:57 2021

@author: grahamseasons
"""

def group_construction(subjects, split_half=True, average=False):
        import itertools
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
            
        return group_container


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
        
copes = [['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_04/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz']]#[['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c1/cope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c2/cope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c3/cope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c4/cope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c5/cope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_c/mapflow/_applywarp_c6/cope7_warp.nii.gz']]
groups = [[['01', '02'], ['03', '04']], [['01', '03'], ['02', '04']], [['01', '04'], ['02', '03']], [['02', '03'], ['01', '04']], [['02', '04'], ['01', '03']], [['03', '04'], ['01', '02']]]#[[['03'], ['01', '02']], [['02'], ['01', '03']], [['01'], ['02', '03']]]
varcopes = [['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v1/varcope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v2/varcope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v3/varcope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v4/varcope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v5/varcope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_01/applywarp_v/mapflow/_applywarp_v6/varcope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v1/varcope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v2/varcope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v3/varcope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v4/varcope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v5/varcope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_02/applywarp_v/mapflow/_applywarp_v6/varcope7_warp.nii.gz'], ['/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v1/varcope2_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v2/varcope3_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v3/varcope4_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v4/varcope5_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v5/varcope6_warp.nii.gz', '/Users/grahamseasons/fMRI/output_comp/working_dir/data_proc/warpflow/_subject_03/applywarp_v/mapflow/_applywarp_v6/varcope7_warp.nii.gz']]

#construction(groups, copes, varcopes)
group_construction(['01','02','03','04'])