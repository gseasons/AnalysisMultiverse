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
from updated.l2_analysis.analyze import level2
from updated.l3_analysis.analyze import level3
from updated.correction.correct import correction
from updated.split_half.split import split
from updated.functions import traverse, metadata, event_grabber, covariate_frame, remove, split_
from copy import deepcopy

class analysis:
    def __init__(self, exp_dir, working_dir, data_dir, out_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        self.data_dir = data_dir
        self.out_dir = out_dir
        
    def construct(self, subjects, sessions, runs, task, pipeline, master, dynamic, split_half):
        fmri = Workflow('fmri')
        base_dir = opj(self.exp_dir, self.working_dir)
        fmri.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['mask', 'task']), name='inputnode')
        inputnode.synchronize = True
        
        bids_dg = Node(BIDSDataGrabber(), name='bids_dg')
        bids_dg.inputs.base_dir = self.data_dir
        if type(subjects) == str:
            subjects = [subjects]
        bids_dg.iterables = ('subject', subjects)
        
        pre = preprocess(task, pipeline, self.out_dir)
        pre = pre.construct(dynamic)
        
        l1 = level1(task, pipeline, self.out_dir, networks=1)
        l1 = l1.construct(dynamic)
        
        l2 = level2(task, pipeline, self.out_dir)
        l2 = l2.construct()
        
        l3 = level3(task, pipeline, self.out_dir)
        l3 = l3.construct()
        l3.inputs.inputnode.subjects = subjects
        
        correct = correction(task, pipeline, self.out_dir)
        correct = correct.construct()
        
        join_sub = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_sub', 
                                   joinsource='bids_dg', joinfield=['copes', 'varcopes', 'bold_warped'])

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
            
            join_scans = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_scans', 
                                  joinsource='ses_run', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask')]),
                          (ses_run, get_files, [('sessions', 'session'),
                                                ('runs', 'run')]),
                          (l1, join_scans, [('outnode.cope', 'copes'),
                                            ('outnode.varcope', 'varcopes'),
                                            ('outnode.bold', 'bold_warped')]),
                          (join_scans, l2, [('copes', 'inputnode.copes'),
                                            ('varcopes', 'inputnode.varcopes')]),
                          (l2, join_sub, [('outnode.copes', 'copes'),
                                          ('outnode.varcopes', 'varcopes')]),
                          (join_scans, join_sub, [('bold_warped', 'bold_warped')]),
                          ])
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
            
            join_scans = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_scans', 
                                  joinsource='ses', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask')]),
                          (ses, get_files, [('sessions', 'session')]),
                          (l1, join_scans, [('outnode.cope', 'copes'),
                                            ('outnode.varcope', 'varcopes'),
                                            ('outnode.bold', 'bold_warped')]),
                          (join_scans, l2, [('copes', 'inputnode.copes'),
                                            ('varcopes', 'inputnode.varcopes')]),
                          (l2, join_sub, [('outnode.copes', 'copes'),
                                          ('outnode.varcopes', 'varcopes')]),
                          (join_scans, join_sub, [('bold_warped', 'bold_warped')]),
                          ])
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
            
            join_scans = JoinNode(IdentityInterface(fields=['copes', 'varcopes', 'bold_warped']), name='join_scans', 
                                  joinsource='run', joinfield=['copes', 'varcopes', 'bold_warped'])
            
            fmri.connect([(inputnode, l2, [('mask', 'inputnode.mask')]),
                          (run, get_files, [('runs', 'run')]),
                          (l1, join_scans, [('outnode.cope', 'copes'),
                                            ('outnode.varcope', 'varcopes'),
                                            ('outnode.bold', 'bold_warped')]),
                          (join_scans, l2, [('copes', 'inputnode.copes'),
                                            ('varcopes', 'inputnode.varcopes')]),
                          (l2, join_sub, [('outnode.copes', 'copes'),
                                          ('outnode.varcopes', 'varcopes')]),
                          (join_scans, join_sub, [('bold_warped', 'bold_warped')]),
                          ])
        else:
            def query(task):
                query = {'bold_files': dict(datatype='func', suffix='bold', task=task, extension=['.nii.gz', '.nii']),
                         'T1w_files': dict(datatype='anat', suffix='T1w', extension=['.nii.gz', '.nii']),
                         }
                return query
        
            get_files = Node(Function(input_names=['task'],
                                      output_names=['query'], function=query), name='get_files')
            
            fmri.connect([(l1, join_sub, [('outnode.cope', 'copes'),
                                          ('outnode.varcope', 'varcopes'),
                                          ('outnode.bold', 'bold_warped')]),
                          ])
        
        meta = Node(Function(input_names=['filename', 'data'], 
                             output_names=['TR'], function=metadata), name='meta')
        meta.inputs.data = self.data_dir
        
        events = Node(Function(input_names=['file', 'data'],
                               output_names=['events'], function=event_grabber), name='events')
        events.inputs.data = self.data_dir
        
        frame = Node(Function(input_names=['data'],
                              output_names=['frame'], function=covariate_frame), name='frame')
        frame.inputs.data = self.data_dir
        
        remove_containers = Node(Function(input_names=['T1w', 'bold'],
                                         output_names=['T1w', 'bold'], function=remove), name='remove_containers')
        
        fmri.connect([(inputnode, pre, [('mask', 'inputnode.mask')]),
                      (inputnode, l1, [('mask', 'inputnode.mask')]),
                      (inputnode, l3, [('mask', 'inputnode.mask')]),
                      (inputnode, get_files, [('task', 'task')]),
                      (get_files, bids_dg, [('query', 'output_query')]),
                      (bids_dg, remove_containers, [('T1w_files', 'T1w'),
                                                    ('bold_files', 'bold')]),
                      (remove_containers, pre, [('T1w', 'inputnode.T1w'),
                                                ('bold', 'inputnode.bold')]),
                      (remove_containers, meta, [('bold', 'filename')]),
                      (remove_containers, events, [('bold', 'file')]),
                      (meta, pre, [('TR', 'inputnode.TR')]),
                      (meta, l1, [('TR', 'inputnode.TR')]),
                      (pre, l1, [('outnode.segmentations', 'inputnode.segmentations'),
                                 ('outnode.warp_file', 'inputnode.warp_file'),
                                 #('outnode.outliers', 'inputnode.outliers'),
                                 ('outnode.brain', 'inputnode.brain'),
                                 ('outnode.brainmask', 'inputnode.brainmask'),
                                 ('outnode.invwarp', 'inputnode.invwarp')]),
                                 #('outnode.mc_par', 'inputnode.mc_par')]),
                      (events, l1, [('events', 'inputnode.event_file')]),
                      (frame, l3, [('frame', 'inputnode.covariates')]),
                      (join_sub, l3, [('copes', 'inputnode.copes'),
                                      ('varcopes', 'inputnode.varcopes')]),
                      (l3, correct, [('outnode.copes', 'inputnode.copes'),
                                     ('outnode.zstats', 'inputnode.zstat'),
                                     ('outnode.mask', 'inputnode.mask')]),
                      ])
        
        traverse(master, fmri, '_full')
        
        if split_half:
            splithalf = split(task, pipeline, self.out_dir)
            splithalf = splithalf.construct()
            splitflow = deepcopy(fmri)
            splitflow.name = 'splitflow'
            if hasattr(splitflow, '_id'):
                splitflow._id = 'splitflow'
                
            divide = Node(Function(input_names=['smoothed', 'unsmoothed', 'half'],
                                   output_names=['smoothed', 'unsmoothed'], function=split_), name='divide')
            divide.iterables = ('half', ['first', 'second'])
            
            join_bold = JoinNode(IdentityInterface(fields='bold'), name='join_bold', 
                                 joinsource='splitflow.l1.outnode', joinfield='bold')
            join_corrected = JoinNode(IdentityInterface(fields=['corrected']), name='join_corrected', 
                                      joinsource='splitflow.correct.outnode', joinfield='corrected')
            
            fmri.connect([(frame, splithalf, [('frame', 'inputnode.covariates')]),
                          (pre, divide, [('outnode.smoothed', 'smoothed'),
                                         ('outnode.unsmoothed', 'unsmoothed')]),
                          (divide, splitflow.get_node('l1'), [('smoothed', 'inputnode.smoothed'),
                                                              ('unsmoothed', 'inputnode.unsmoothed')]),
                          (splitflow.get_node('l1'), join_bold, [('outnode.bold', 'bold')]),
                          (splitflow.get_node('l3'), splithalf, [('outnode.mask', 'mask')]),
                          (splitflow.get_node('corrected'), join_corrected, [('outnode.corrected_files', 'corrected')]),
                          (join_bold, splithalf, [('bold', 'bold')]),
                          (join_corrected, splithalf, [('corrected', 'corrected')]),
                          ])
            
        fmri.connect([(pre, l1, [('outnode.smoothed', 'inputnode.smoothed'),
                                 ('outnode.unsmoothed', 'inputnode.unsmoothed')]),
                      ])
        
        return fmri