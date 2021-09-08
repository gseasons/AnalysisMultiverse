#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 10:41:22 2021

@author: grahamseasons
"""
from workflows import *
from iterface import IdentityIterface
#NEED TO ACCOUNT FOR VARIABLE RUNS AND SESSIONS ACROSS TASKS
#SEND DESIRED INTERMEDIATES AND OUTPUTS TO DATASINK

class analysis:
    def __init__(self, exp_dir, working_dir, data_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        self.data_dir = data_dir
        
    def construct(self, subjects, tasks, sessions, runs, num_pipelines):
        fmri = Workflow('fmri')
        base_dir = opj(self.exp_dir, self.working_dir)
        fmri.base_dir = base_dir
        
        inputnode = Node(IdentityIterface(fields=['pipeline_names',
                                                  'frac_mask',
                                                  'discard',
                                                  'dof_mc',
                                                  'fwhm',
                                                  'cost_mc',
                                                  'bet_frac',
                                                  'robust',
                                                  'wm_thresh',
                                                  'dof_f',
                                                  'bbr_type',
                                                  'interp',
                                                  'cost',
                                                  'bins',
                                                  'iso',
                                                  'bbr',
                                                  'susan',
                                                  'HP',
                                                  'film_thresh',
                                                  'serial_cor',
                                                  'base_switch',
                                                  'gamma',
                                                  'dict_opt',
                                                  'base_val',
                                                  'mask',
                                                  'mode_l2',
                                                  'mode_l3',
                                                  'method',
                                                  'p',
                                                  'connectivity',
                                                  'z_thresh',
                                                  'index']), name='inputnode')
        
        
        outnode = Node(IdentityInterface(fields=['TBD']), name='outnode')
        
        iternode = Node(IdentityInterface(fields=['pipeline']), name='iternode')
        iternode.iterables = [('pipeline', range(num_pipelines))]
        
        task = Node(IdentityInterface(fields=['task']), name='task')
        task.iterables = [('task', tasks)]
        
        bids_dg = Node(BIDSDataGrabber(), name='bids_dg')
        bids_dg.inputs.base_dir = self.data_dir
        bids_dg.iterables = ('subject', subjects)
        
        pre = preprocess(base_dir).construct()
        
        l1 = level1(base_dir).construct()
        
        norm = spatial_normalization(base_dir).construct()
        
        join_sub = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_sub', 
                                   joinsource='bids_dg', joinfield=['copes', 'varcopes', 'bold_warped'])
        
        l2_split = level2(base_dir).construct(len(subjects))
        l2_split.inputs.inputnode.subjects = subjects
        l2_split.inputs.inputnode.split_half = True
        split = split_half(base_dir).construct()
        
        l3 = level3(base_dir).construct()
        l3.inputs.inputnode.subjects = subjects
        
        correct = correction(base_dir).construct()
        
        def remove_container(file):
            return file[0]
        
        #Probably replace session and run checks with MapNodes
        if sessions and runs:
            ses_run = Node(IdentityInterface(fields=['sessions', 'runs']), name='ses_run')
            ses_run.iterables = [('sessions', sessions),
                                 ('runs', runs)]
            
            def query(session, run, task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, session=session, run=run),
                         'T1w_files': dict(datatype='anat', suffix='T1w', session=session),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['session', 'run', 'task'],
                                      output_names=['query'], function=query), name='get_files')
            
            join_sesrun = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_sesrun', 
                                   joinsource='ses_run', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            l2 = level2(base_dir).construct(len(subjects), '_fe')
            l2.inputs.inputnode.subjects = subjects
            l2.inputs.inputnode.split_half = False
            
            fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask'),
                                           ('mode_l2', 'inputnode.mode')]),
                          (ses_run, get_files, [('sessions', 'session'),
                                                ('runs', 'run')]),
                          (norm, join_sesrun, [('outnode.cope', 'copes'),
                                               ('outnode.varcope', 'varcopes'),
                                               ('outnode.bold', 'bold_warped')]),
                          (join_sesrun, l2, [('copes', 'inputnode.copes'),
                                             ('varcopes', 'inputnode.varcopes')]),
                          (l2, join_sub, [(('outnode.copes', remove_container), 'copes'),
                                          (('outnode.varcopes', remove_container), 'varcopes')]),
                          (join_sesrun, join_sub, [('bold_warped', 'bold_warped')]),
                          ])
        elif sessions:
            ses = Node(IdentityInterface(fields=['sessions']), name='ses')
            ses.iterables = [('sessions', sessions)]
            
            def query(session, task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, session=session),
                         'T1w_files': dict(datatype='anat', suffix='T1w', session=session),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['session', 'task'],
                                      output_names=['query'], function=query), name='get_files')
            
            join_ses = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_ses', 
                                   joinsource='ses', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            l2 = level2(base_dir).construct(len(subjects), '_fe')
            l2.inputs.inputnode.subjects = subjects
            l2.inputs.inputnode.split_half = False
            
            fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask'),
                                           ('mode_l2', 'inputnode.mode')]),
                          (ses, get_files, [('sessions', 'session')]),
                          (norm, join_ses, [('outnode.cope', 'copes'),
                                            ('outnode.varcope', 'varcopes'),
                                            ('outnode.bold', 'bold_warped')]),
                          (join_ses, l2, [('copes', 'inputnode.copes'),
                                          ('varcopes', 'inputnode.varcopes')]),
                          (l2, join_sub, [(('outnode.copes', remove_container), 'copes'),
                                          (('outnode.varcopes', remove_container), 'varcopes')]),
                          (join_ses, join_sub, [('bold_warped', 'bold_warped')]),
                          ])
        elif runs:
            run = Node(IdentityInterface(fields=['runs']), name='run')
            run.iterables = [('runs', runs)]
            
            def query(run, task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, run=run),
                         'T1w_files': dict(datatype='anat', suffix='T1w'),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['run', 'task'],
                                      output_names=['query'], function=query), name='get_files')
            
            join_run = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_run', 
                                   joinsource='run', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            l2 = level2(base_dir).construct(len(subjects), '_fe')
            l2.inputs.inputnode.subjects = subjects
            l2.inputs.inputnode.split_half = False
            
            fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask'),
                                           ('mode_l2', 'inputnode.mode')]),
                          (run, get_files, [('runs', 'run')]),
                          (norm, join_run, [('outnode.cope', 'copes'),
                                            ('outnode.varcope', 'varcopes'),
                                            ('outnode.bold', 'bold_warped')]),
                          (join_run, l2, [('copes', 'inputnode.copes'),
                                          ('varcopes', 'inputnode.varcopes')]),
                          (l2, join_sub, [(('outnode.copes', remove_container), 'copes'),
                                          (('outnode.varcopes', remove_container), 'varcopes')]),
                          (join_run, join_sub, [('bold_warped', 'bold_warped')]),
                          ])
        else:
            def query(task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task),
                         'T1w_files': dict(datatype='anat', suffix='T1w'),
                         }
                return query
                
            
            get_files = Node(Function(input_names=['task'],
                                      output_names=['query'], function=query), name='get_files')
            
            fmri.connect([(norm, join_sub, [('outnode.cope', 'copes'),
                                            ('outnode.varcope', 'varcopes'),
                                            ('outnode.bold', 'bold_warped')]),
                          ])
            
            
        
        #NOTE: NEED TO DEAL WITH CASE WHERE SESSIONS AND RUNS ARE DEPENDENT ON TASK
        #TR, event_file, mc_par, outliers, smoothed -> l1
        #covariate_frame
        
        fmri.connect([(iternode, inputnode, [('pipeline', 'index')]),
                      (inputnode, pre, [('frac_mask', 'inputnode.frac_mask'),
                                        ('discard', 'inputnode.discard'),
                                        ('dof_mc', 'inputnode.dof_mc'),
                                        ('fwhm', 'inputnode.fwhm'),
                                        ('cost_mc', 'inputnode.cost_mc'),
                                        ('bet_frac', 'inputnode.bet_frac'),
                                        ('robust', 'inputnode.robust'),
                                        ('wm_thresh', 'inputnode.wm_thresh'),
                                        ('dof_f', 'inputnode.dof_f'),
                                        ('bbr_type', 'inputnode.bbr_type'),
                                        ('interp', 'inputnode.interp'),
                                        ('cost', 'inputnode.cost'),
                                        ('bins', 'inputnode.bins'),
                                        ('iso', 'inputnode.iso'),
                                        ('bbr', 'inputnode.bbr'),
                                        ('susan', 'inputnode.susan')]),
                      (inputnode, l1, [('discard', 'inputnode.discard'),
                                       ('HP', 'inputnode.HP'),
                                       ('film_thresh', 'inputnode.thresh'),
                                       ('serial_cor', 'inputnode.serial_cor'),
                                       ('base_switch', 'inputnode.base_switch'),
                                       ('gamma', 'inputnode.gamma'),
                                       ('dict_opt', 'inputnode.dict_opt'),
                                       ('base_val', 'inputnode.base_val')]),
                      (inputnode, norm, [('mask', 'inputnode.ref_file')]),
                      (inputnode, l2_split, [('mode_l2', 'inputnode.mode'),
                                             ('mask', 'inputnode.mask')]),
                      (inputnode, split, [('mask', 'inputnode.mask')]),
                      (inputnode, l3, [('mask', 'inputnode.mask'),
                                       ('mode_l3', 'inputnode.mode')]),
                      (inputnode, correct, [('method', 'inputnode.method'),
                                            ('p', 'inputnode.p'),
                                            ('mask', 'inputnode.mask'),
                                            ('connectivity', 'inputnode.connectivity'),
                                            ('z_thresh', 'inputnode.z_thresh')]),
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
            event_file = layout.get(task=task, extension='.tsv')
            
            if len(event_file) > 1:
                sub = re.search('/sub-([0-9]+)/', file).group(1)
                ses = re.search('_ses-([A-Za-z]+)_task', file)
                run = re.search('run-([0-9]+)', file)
                
                if ses and run:
                    event_file = layout.get(task=task, session=ses.group(1), run=run.group(1), subject=sub, extension='.tsv')
                elif ses:
                    event_file = layout.get(task=task, session=ses.group(1), subject=sub, extension='.tsv')
                elif run:
                    event_file = layout.get(task=task, run=run.group(1), subject=sub, extension='.tsv')
                
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
        
        join_task = JoinNode(IdentityInterface(fields=['score']), name='join_task', joinsource='task', joinfield=['score'])
        
        
        fmri.connect([(task, get_files, [('task', 'task')]),
                      (get_files, bids_dg, [('query', 'output_query')]),
                      #(bids_dg, find_bold, [('bold_files', 'files'),
                       #                     ('task', 'task')]),
                      #(find_bold, pre, [('file', 'inputnode.bold')]),
                      (bids_dg, pre, [(('T1w_files', remove_container), 'inputnode.T1w'),
                                      (('bold_files', remove_container), 'inputnode.bold')]),
                      (bids_dg, meta, [(('bold_files', remove_container), 'filename')]),
                      (bids_dg, events, [(('bold_files', remove_container), 'file')]),
                      (meta, pre, [('TR', 'inputnode.TR')]),
                      (meta, l1, [('TR', 'inputnode.TR')]),
                      (events, l1, [('events', 'inputnode.event_file')]),
                      (pre, l1, [('outnode.smoothed', 'inputnode.smoothed'),
                                 ('outnode.outliers', 'inputnode.outliers'),
                                 ('outnode.mc_par', 'inputnode.mc_par')]),
                      (l1, norm, [('outnode.feat_dir', 'inputnode.feat_dir')]),
                      (pre, norm, [('outnode.brain', 'inputnode.brain')]),
                      (join_sub, l2_split, [('copes', 'inputnode.copes'),
                                            ('varcopes', 'inputnode.varcopes')]),
                      (l2_split, split, [('outnode.groups', 'inputnode.groups'),
                                          ('outnode.zstats', 'inputnode.zstats')]),
                      (join_sub, split, [('bold_warped', 'inputnode.preproc_bold')]),
                      (frame, split, [('frame', 'inputnode.covariate_frame')]),
                      (frame, l3, [('frame', 'inputnode.covariates')]),
                      (join_sub, l3, [('copes', 'inputnode.copes'),
                                      ('varcopes', 'inputnode.varcopes')]),
                      (l3, correct, [('outnode.copes', 'inputnode.copes'),
                                        ('outnode.zstats', 'inputnode.zstat')]),
                      (split, join_task, [('outnode.score', 'score')]),
                      ])
        
        return fmri
    
#A = analysis('/Users/grahamseasons/fMRI/output_comp', 'working_dir', "/Users/grahamseasons/fMRI/test_data/ds000114/").construct(['01','02'], ['fingerfootlips'], ['test'], [], 1)
#B = 3
            

            
        
        