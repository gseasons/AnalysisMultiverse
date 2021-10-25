#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 13:31:06 2021

@author: grahamseasons
"""
from nipype import Workflow, Node, JoinNode, IdentityInterface, Function
from nipype.interfaces.io import BIDSDataGrabber
from os.path import join as opj
from updated.preprocessing.preprocess import preprocess
from updated.l1_analysis.analyze import level1
from updated.functions import traverse

class analysis:
    def __init__(self, exp_dir, working_dir, data_dir, out_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        self.data_dir = data_dir
        self.out_dir = out_dir
        
    def construct(self, subjects, sessions, runs, task, pipeline, master, dynamic):
        fmri = Workflow('fmri')
        base_dir = opj(self.exp_dir, self.working_dir)
        fmri.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['mask', 'task', 'max_groups']), name='inputnode')
        inputnode.synchronize = True
        
        
        outnode = Node(IdentityInterface(fields=['TBD']), name='outnode')
        
        #iternode = Node(IdentityInterface(fields=['pipeline']), name='iternode')
        #iternode.iterables = [('pipeline', range(num_pipelines))]
        
        bids_dg = Node(BIDSDataGrabber(), name='bids_dg')
        bids_dg.inputs.base_dir = self.data_dir
        if type(subjects) == str:
            subjects = [subjects]
        bids_dg.iterables = ('subject', subjects)
        
        pre = preprocess(task, pipeline, self.out_dir)
        pre = pre.construct(dynamic)
        
        l1 = level1(task, pipeline, self.out_dir, networks=1)
        l1 = l1.construct(dynamic)
        #ADD ALL OTHER NODES HERE TOO
        #MAKE CONNECTIONS ONCE VERIFY TRAVERSE WORKING AS EXPECTED FOR L1
        fmri.add_nodes([pre, l1])
        traverse(master, fmri, '_full')
        
        #l1 = level1(self.out_dir).construct(l1_v, spatial_norm)
        
        #norm = spatial_normalization(base_dir)
        #norm_get = norm.construct()
        #norm_app = norm.construct()
        
        join_sub = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_sub', 
                                   joinsource='bids_dg', joinfield=['copes', 'varcopes', 'bold_warped'])
        
        #l2_split = level2(self.out_dir).construct(group_num, l2_v)
        #l2_split.inputs.inputnode.subjects = subjects
        #l2_split.inputs.inputnode.split_half = True
        #split = split_half(self.out_dir).construct(pipeline_ind)
        
        #l3 = level3(self.out_dir).construct(l3_v)
        #l3.inputs.inputnode.subjects = subjects
        
        #correct = correction(self.out_dir).construct(corr)
        def remove_container(file):
            return file[0]
        
        def ident(file):
            return file
        
        def fix_data_struct(file):
            return file
        
        #Probably replace session and run checks with MapNodes
        if sessions and runs:
            ses_run = Node(IdentityInterface(fields=['sessions', 'runs']), name='ses_run')
            ses_run.iterables = [('sessions', sessions),
                                 ('runs', runs)]
            
            def query(session, run, task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, session=session, run=run, extension=['.nii.gz', '.nii']),
                         'T1w_files': dict(datatype='anat', suffix='T1w', session=session, extension=['.nii.gz', '.nii']),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['session', 'run', 'task'],
                                      output_names=['query'], function=query), name='get_files')
            
            join_sesrun = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_sesrun', 
                                   joinsource='ses_run', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            #l2 = level2(base_dir).construct(group_num, l2_v, '_fe')
            #l2.inputs.inputnode.subjects = subjects
            #l2.inputs.inputnode.split_half = False
            
# =============================================================================
#             fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask'),
#                                            ('max_groups', 'inputnode.max')]),
#                           (ses_run, get_files, [('sessions', 'session'),
#                                                 ('runs', 'run')]),
#                           (l1, join_sesrun, [('outnode.cope', 'copes'),
#                                                ('outnode.varcope', 'varcopes'),
#                                                ('outnode.bold', 'bold_warped')]),
#                           (join_sesrun, l2, [('copes', 'inputnode.copes'),
#                                              ('varcopes', 'inputnode.varcopes')]),
#                           (l2, join_sub, [(('outnode.copes', ident), 'copes'),
#                                           (('outnode.varcopes', ident), 'varcopes')]),
#                           (join_sesrun, join_sub, [('bold_warped', 'bold_warped')]),
#                           ])
# =============================================================================
        elif sessions:
            ses = Node(IdentityInterface(fields=['sessions']), name='ses')
            ses.iterables = [('sessions', sessions)]
            
            def query(session, task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, session=session, extension=['.nii.gz', '.nii']),
                         'T1w_files': dict(datatype='anat', suffix='T1w', session=session, extension=['.nii.gz', '.nii']),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['session', 'task'],
                                      output_names=['query'], function=query), name='get_files')
            #PERHAPS CHANGE SO CAN RUN ON ONLY ONE SUBJECT OR SESSION
            join_ses = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_ses', 
                                   joinsource='ses', joinfield=['copes', 'varcopes', 'bold_warped'])
            
# =============================================================================
#             l2 = level2(base_dir).construct(group_num, l2_v, '_fe')
#             l2.inputs.inputnode.subjects = subjects
#             l2.inputs.inputnode.split_half = False
#             
#             fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask'),
#                                            ('max_groups', 'inputnode.max')]),
#                           (ses, get_files, [('sessions', 'session')]),
#                           (l1, join_ses, [('outnode.cope', 'copes'),
#                                             ('outnode.varcope', 'varcopes'),
#                                             ('outnode.bold', 'bold_warped')]),
#                           (join_ses, l2, [('copes', 'inputnode.copes'),
#                                           ('varcopes', 'inputnode.varcopes')]),
#                           (l2, join_sub, [(('outnode.copes', ident), 'copes'),
#                                           (('outnode.varcopes', ident), 'varcopes')]),
#                           (join_ses, join_sub, [('bold_warped', 'bold_warped')]),
#                           ])
# =============================================================================
        elif runs:
            run = Node(IdentityInterface(fields=['runs']), name='run')
            run.iterables = [('runs', runs)]
            
            def query(run, task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, run=run, extension=['.nii.gz', '.nii']),
                         'T1w_files': dict(datatype='anat', suffix='T1w', extension=['.nii.gz', '.nii']),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['run', 'task'],
                                      output_names=['query'], function=query), name='get_files')
            
            join_run = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_run', 
                                   joinsource='run', joinfield=['copes', 'varcopes', 'bold_warped'])
            
# =============================================================================
#             l2 = level2(base_dir).construct(group_num, l2_v, '_fe')
#             l2.inputs.inputnode.subjects = subjects
#             l2.inputs.inputnode.split_half = False
#             
#             fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask'),
#                                            ('max_groups', 'inputnode.max')]),
#                           (run, get_files, [('runs', 'run')]),
#                           (l1, join_run, [('outnode.cope', 'copes'),
#                                             ('outnode.varcope', 'varcopes'),
#                                             ('outnode.bold', 'bold_warped')]),
#                           (join_run, l2, [('copes', 'inputnode.copes'),
#                                           ('varcopes', 'inputnode.varcopes')]),
#                           (l2, join_sub, [(('outnode.copes', ident), 'copes'),
#                                           (('outnode.varcopes', ident), 'varcopes')]),
#                           (join_run, join_sub, [('bold_warped', 'bold_warped')]),
#                           ])
# =============================================================================
        else:
            def query(task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, extension=['.nii.gz', '.nii']),
                         'T1w_files': dict(datatype='anat', suffix='T1w', extension=['.nii.gz', '.nii']),
                         }
                return query
                
            def fix_data_struct(file):
                files = []
                if len(file) > 1:
                    for pipeline in file:
                        files.append([pipeline])
                else:
                    files = file
                return files
        
            get_files = Node(Function(input_names=['task'],
                                      output_names=['query'], function=query), name='get_files')
            
# =============================================================================
#             fmri.connect([(l1, join_sub, [(('outnode.cope', fix_data_struct), 'copes'),
#                                             (('outnode.varcope', fix_data_struct), 'varcopes'),
#                                             (('outnode.bold', fix_data_struct), 'bold_warped')]),
#                           ])
# =============================================================================
            
            
        
        #NOTE: NEED TO DEAL WITH CASE WHERE SESSIONS AND RUNS ARE DEPENDENT ON TASK
        #TR, event_file, mc_par, outliers, smoothed -> l1
        #covariate_frame
        
        fmri.connect([#(iternode, inputnode, [('pipeline', 'index')]),
                      (inputnode, pre, [#('frac_mask', 'inputnode.frac_mask'),
                                        #('discard', 'inputnode.discard'),
                                        #('dof_mc', 'inputnode.dof_mc'),
                                        #('fwhm', 'inputnode.fwhm'),
                                        #('cost_mc', 'inputnode.cost_mc'),
                                        #('bet_frac', 'inputnode.bet_frac'),
                                        #('robust', 'inputnode.robust'),
                                        #('wm_thresh', 'inputnode.wm_thresh'),
                                        #('dof_f', 'inputnode.dof_f'),
                                        #('bbr_type', 'inputnode.bbr_type'),
                                        #('interp', 'inputnode.interp'),
                                        #('cost', 'inputnode.cost'),
                                        #('bins', 'inputnode.bins'),
                                        #('iso', 'inputnode.iso'),
                                        #('bbr', 'inputnode.bbr'),
                                        #('susan', 'inputnode.susan'),
                                        #('warp_post_feat', 'inputnode.warplater'),
                                        ('mask', 'inputnode.mask')]),
                      #(norm_get, pre, [('outnode.warp', 'inputnode.warp_file')]),
                      
                      #(pre, l1, [('outnode.invwarp', 'inputnode.invwarp'),
                      #           ('outnode.warp_file', 'inputnode.warp')]),
# =============================================================================
                     (inputnode, l1, [#('task', 'inputnode.task'),
                                      ('mask', 'inputnode.mask')]),
# =============================================================================
                      
                      #(inputnode, l1, [#('discard', 'inputnode.discard'),
                                       #('HP', 'inputnode.HP'),
                                       #('film_thresh', 'inputnode.thresh'),
                                       #('serial_cor', 'inputnode.serial_cor'),
                                       #('base_switch', 'inputnode.base_switch'),
                                       #('gamma', 'inputnode.gamma'),
                                       #('dict_opt', 'inputnode.dict_opt'),
                                       #('base_val', 'inputnode.base_val'),
                                       #('resting', 'inputnode.resting'),
                                       #('mask', 'inputnode.mask'),
                        #               ('warp_post_feat', 'inputnode.warp_post_feat')]),
                      #(inputnode, norm_get, [('mask', 'inputnode.ref_file')]),
                      #(inputnode, norm_app, [('mask', 'inputnode.ref_file')]),
                      #(inputnode, norm_app, [('warp_post_feat', 'inputnode.needwarp')]),
                      #(norm_get, norm_app, [('outnode.warp', 'inputnode.warp_file')]),
# =============================================================================
#                       (inputnode, l2_split, [#('mode_l2', 'inputnode.mode'),
#                                              ('mask', 'inputnode.mask'),
#                                              ('max_groups', 'inputnode.max')]),
#                       (inputnode, split, [('mask', 'inputnode.mask')]),
#                       (inputnode, l3, [('mask', 'inputnode.mask')]),
# =============================================================================
# =============================================================================
#                       (inputnode, correct, [#('method', 'inputnode.method'),
#                                             #('p', 'inputnode.p'),
#                                             ('mask', 'inputnode.mask')])
# =============================================================================
                                            #('connectivity', 'inputnode.connectivity'),
                                            #('z_thresh', 'inputnode.z_thresh')]),
                      ])
        
        
        def metadata(filename, data):
            from bids.layout import BIDSLayout
            layout = BIDSLayout(data)#"/Users/grahamseasons/fMRI/test_data/ds000114/")
            #meta_path = layout.get(task=task, session=ses, datatype='func', suffix='bold', return_type='filename')[0]
            meta_data = layout.get_metadata(filename)
            TR = meta_data['RepetitionTime']
            return TR
        
        meta = Node(Function(input_names=['filename', 'data'], 
                             output_names=['TR'], function=metadata), name='meta')
        meta.inputs.data = self.data_dir
        
        def event_grabber(file, data):
            import re
            from bids.layout import BIDSLayout
            layout = BIDSLayout(data)
            
            task = re.search('task-([0-9A-Za-z]+)_bold', file).group(1)
            
            if 'rest' in task:
                return ''
            
            event_file = layout.get(task=task, extension='.tsv')
            
            if len(event_file) > 1:
                sub = re.search('/sub-([0-9A-Za-z]+)/', file).group(1)
                ses = re.search('_ses-([A-Za-z]+)_task', file)
                run = re.search('run-([0-9]+)', file)
                
                if ses and run:
                    event_file = layout.get(task=task, session=ses.group(1), run=run.group(1), subject=sub, extension='.tsv')
                elif ses:
                    event_file = layout.get(task=task, session=ses.group(1), subject=sub, extension='.tsv')
                elif run:
                    event_file = layout.get(task=task, run=run.group(1), subject=sub, extension='.tsv')
            elif not len(event_file):
                event_file = ['']
                
            return event_file[0]
        
        events = Node(Function(input_names=['file', 'data'],
                               output_names=['events'], function=event_grabber), name='events')
        events.inputs.data = self.data_dir
        
        def covariate_frame(data):
            from bids.layout import BIDSLayout
            import pandas as pd
            layout = BIDSLayout(data)
            file = layout.get(return_type='filename', extension='.tsv', suffix='participants')
            #frame = pd.read_table(file[0])
            
            return file[0]
        
        frame = Node(Function(input_names=['data'],
                              output_names=['frame'], function=covariate_frame), name='frame')
        frame.inputs.data = self.data_dir
        
        #pre outnode 'smoothed', 'outliers', 'plots', 'mc_par',
                                              #   'warped_mean', 'warped', 'brain', 'reg_out_mat'
        
        #l1 outnode feat_dir, contrast_names
        
        #l2 outnode groups, copes, varcopes, flameo_stats, zstats
        
        #split outnode score, R, R_lst, P, P_lst
        
        #l3 outnode copes var_copes zstats flameo_stats
        
        def get_bold(files, task):
            for file in files:
                if file == task:
                    return file
                
        find_bold = Node(Function(input_names=['files', 'task'],
                                  output_names=['file'], function=get_bold), name='find_bold')
        
        def remove(T1w, bold):
            return T1w[0], bold[0]
        
        #join_task = JoinNode(IdentityInterface(fields=['score']), name='join_task', joinsource='task', joinfield=['score'])
        remove_containers = Node(Function(input_names=['T1w', 'bold'],
                                         output_names=['T1w', 'bold'], function=remove), name='remove_containers')
        
        fmri.connect([(inputnode, get_files, [('task', 'task')]),
                      (get_files, bids_dg, [('query', 'output_query')]),
                      #(bids_dg, find_bold, [('bold_files', 'files'),
                       #                     ('task', 'task')]),
                      #(find_bold, pre, [('file', 'inputnode.bold')]),
                      (bids_dg, remove_containers, [('T1w_files', 'T1w'),
                                                    ('bold_files', 'bold')]),
                      (remove_containers, pre, [('T1w', 'inputnode.T1w'),
                                                ('bold', 'inputnode.bold')]),
                      (remove_containers, meta, [('bold', 'filename')]),
                      (remove_containers, events, [('bold', 'file')]),
                      (meta, pre, [('TR', 'inputnode.TR')]),
                      (meta, l1, [('TR', 'inputnode.TR')]),
                      (pre, l1, [('outnode.smoothed', 'inputnode.smoothed'),
                                 ('outnode.unsmoothed', 'inputnode.unsmoothed'),
                                 ('outnode.segmentations', 'inputnode.segmentations'),
                                 ('outnode.warp_file', 'inputnode.warp_file'),
                                 ('outnode.outliers', 'inputnode.outliers'),
                                 ('outnode.brain', 'inputnode.brain'),
                                 ('outnode.brainmask', 'inputnode.brainmask'),
                                 ('outnode.invwarp', 'inputnode.invwarp')]),
# =============================================================================
#                       (meta, l1, [('TR', 'inputnode.TR')]),
                      (events, l1, [('events', 'inputnode.event_file')]),
#                       (pre, l1, [('outnode.smoothed', 'inputnode.smoothed'),
#                                  ('outnode.outliers', 'inputnode.outliers'),
#                                  ('outnode.mc_par', 'inputnode.mc_par'),
#                                  ('outnode.brain', 'inputnode.brain'),
#                                  ('outnode.warp_file', 'inputnode.warp'),
#                                  ('outnode.segmentations', 'inputnode.segmentations')]),
# =============================================================================
                      #(l1, norm_app, [('outnode.feat_dir', 'inputnode.feat_dir')]),
                      #(pre, norm_get, [('outnode.brain', 'inputnode.brain')]),
# =============================================================================
#                       (join_sub, l2_split, [('copes', 'inputnode.copes'),
#                                             ('varcopes', 'inputnode.varcopes')]),
#                       (l2_split, split, [('outnode.groups', 'inputnode.groups'),
#                                           ('outnode.zstats', 'inputnode.zstats')]),
#                       (join_sub, split, [('bold_warped', 'inputnode.preproc_bold')]),
#                       (frame, split, [('frame', 'inputnode.covariate_frame')]),
#                       (frame, l3, [('frame', 'inputnode.covariates')]),
#                       (join_sub, l3, [('copes', 'inputnode.copes'),
#                                       ('varcopes', 'inputnode.varcopes')]),
#                       (l3, correct, [('outnode.copes', 'inputnode.copes'),
#                                         ('outnode.zstats', 'inputnode.zstat')]),
# =============================================================================
                      #(split, join_task, [('outnode.score', 'score')]),
                      ])
        
        return fmri