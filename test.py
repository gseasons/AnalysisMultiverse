#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:58:44 2021

@author: grahamseasons
"""

def write_out(base_dir, pipeline_st, task, cope, varcope, bold, feat_dir, seed):
    outputs = list(locals().keys())[3:]
    from nipype.interfaces import DataSink
    from nipype import Node
    import re, os
    sink = Node(DataSink(base_directory=base_dir, parameterization=True), name='sink')
    folder_name = task
    test_str = vars()[outputs[0]]
    if isinstance(test_str, list):
        test_str = test_str[0]
    session = re.search('/(_sessions_[A-Za-z0-9]+)/', test_str)
    run = re.search('/(_runs_[0-9]+)/', test_str)
    subject = re.search('/(_subject_[0-9A-Za-z]+/)', test_str).group(1)

    if session:
        folder_name += '/' + session.group(1) + '/'
    if run:
        folder_name += '/' + run.group(1) + '/'
    if folder_name == task:
        folder_name += '/'

    #folder_name += wf_name + '/' + subject

    for i, out in enumerate(outputs):
        if out:

            if isinstance(vars()[out], str) and os.path.isdir(vars()[out]):
                        for file in ['design.con', 'design.mat']:
                            file = vars()[out] + '/' + file
                            if os.path.isdir(file) or os.path.isfile(file):
                                setattr(sink.inputs, 'pipelines/' + task + '.@' + str(i) + file, file)

            else:
                setattr(sink.inputs, 'pipelines/' + task + '.@' + str(i), vars()[out])

    sink.run()

base_dir = '/Volumes/NewVolume/sup_pre_l1_v2/processed'
bold = ['/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/level1/_subject_002S6009/_i_0/_i_0/_i_2/_network_0/applywarpbold/filtered_func_data_warp.nii.gz']
cope = ['/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/level1/_subject_002S6009/_i_0/_i_0/_i_2/_network_0/applywarpcopes/mapflow/_applywarpcopes0/cope1_warp.nii.gz']
feat_dir = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/level1/_subject_002S6009/_i_0/_i_0/_i_2/_network_0/feat/run0.feat'
pipeline_st = 0
seed = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/level1/_subject_002S6009/_i_0/_i_0/_i_2/_network_0/Finfo/sub-002S6009_task-rest_bold_roi_mcf_regressed_maths_filt_smooth_reho_thresh_maths.nii.gz'
task = 'rest'
varcope = ['/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/level1/_subject_002S6009/_i_0/_i_0/_i_2/_network_0/applywarpvarcopes/mapflow/_applywarpvarcopes0/varcope1_warp.nii.gz']
write_out(base_dir, pipeline_st, task, cope, varcope, bold, feat_dir, seed)