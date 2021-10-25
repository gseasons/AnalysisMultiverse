#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:36:22 2021

@author: grahamseasons
"""
from updated.normalization.spatial_normalization import spatial_normalization
from nipype import Workflow, Node, Function, IdentityInterface
from versatile import Level1DesignVersatile
from nipype.interfaces.fsl import FEAT
import os

class level1(spatial_normalization):
    def __init__(self, task, pipeline, base_dir, networks):
        super().__init__()
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        self.networks = networks

    def construct(self, func_dic):
        from updated.l1_analysis.functions import get_sink
        level1 = Workflow('level1')
        level1.base_dir = os.getcwd()
        inputnode = Node(IdentityInterface(fields=['smoothed', 'unsmoothed', 'brainmask', 'segmentations', 'warp_file', 'outliers', 'brain', 'invwarp',
                                                   'event_file', 'TR', 'mask']), name='inputnode')
        
        intermediates = ['cope', 'varcope', 'bold', 'feat_dir', 'seed']
        files = ['design.con', 'design.mat']
        
        self.l1(level1, func_dic)
        self.apply_warps(level1)
        
        level1.connect([(inputnode, level1.get_node('Finfo'), [('mask', 'mask'),
                                                               ('event_file', 'event_file'),
                                                               ('TR', 'TR'),
                                                               ('smoothed', 'smoothed'),
                                                               ('unsmoothed', 'unsmoothed'),
                                                               ('brainmask', 'brainmask'),
                                                               ('brain', 'brain'),
                                                               ('segmentations', 'segmentations'),
                                                               ('invwarp', 'invwarp'),
                                                               ('outliers', 'outliers'),
                                                               ]),
                        (inputnode, level1.get_node('l1d'), [('TR', 'interscan_interval')]),
                        (inputnode, level1.get_node('correction'), [('TR', 'TR')]),
                        (inputnode, level1.get_node('applywarpcopes'), [('warp_file', 'field_file')]),
                        (inputnode, level1.get_node('applywarpvarcopes'), [('warp_file', 'field_file')]),
                        (inputnode, level1.get_node('applywarpbold'), [('warp_file', 'field_file')]),
                        (inputnode, level1.get_node('applywarpcopes'), [('mask', 'ref_file')]),
                        (inputnode, level1.get_node('applywarpvarcopes'), [('mask', 'ref_file')]),
                        (inputnode, level1.get_node('applywarpbold'), [('mask', 'ref_file')]),
                        (level1.get_node('feat'), level1.get_node('selectfiles'), [('feat_dir', 'base_directory')]),
                        ])
        
        outnode = Node(IdentityInterface(fields=['feat_dir', 'cope', 'varcope', 'bold', 'ev_files']), name='outnode')
        
        level1.connect([(level1.get_node('feat'), outnode, [('feat_dir', 'feat_dir')]),
                        (level1.get_node('ret'), outnode, [('cope', 'cope'),
                                                           ('varcope', 'varcope'),
                                                           ('bold', 'bold')]),
                        (level1.get_node('l1d'), outnode, [('ev_files', 'ev_files')]),
                        ])
        
        write = Node(Function(input_names=['base_dir', 'pipeline_st', 'task']+intermediates), name='write')
        write.inputs.function_str = get_sink(intermediates, files)
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        level1.connect([(level1.get_node('feat'), write, [('feat_dir', 'feat_dir')]),
                        (level1.get_node('ret'), write, [('cope', 'cope'),
                                                         ('varcope', 'varcope'),
                                                         ('bold', 'bold')]),
                        ])
        
        if 'rest' in self.task:
            level1.connect([(level1.get_node('Finfo'), write, [('seed', 'seed')])])
        
        return level1
        
        
    def l1(self, flow, func_dic):
        from updated.l1_analysis.functions import function_str, correct_task_timing, contrasts
        networks = Node(IdentityInterface(fields=['network']), name='networks')
        networks.iterables = ('network', range(self.networks))
        
        func_str, input_names = function_str('info', func_dic)
        if 'rest' in self.task:
            outnames = ['session_info', 'seed']
        else:
            outnames = ['session_info']
        Finfo = Node(Function(input_names=input_names, output_names=outnames), name='Finfo')
        Finfo.inputs.function_str = func_str
        Finfo.inputs.task = self.task
        
        correction = Node(Function(input_names=['session_info', 'TR', 'discard'], 
                                   output_names=['session_info'], 
                                   function=correct_task_timing), name='correction')
        
        contrasts = Node(Function(input_names=['session_info'], 
                                  output_names=['identities', 'contrasts'],
                                  function=contrasts), name='contrasts')
        
        l1d = Node(Level1DesignVersatile(), name='l1d')
        feat = Node(FEAT(), name='feat')
        
        flow.connect([(networks, Finfo, [('network', 'network')]),
                      (Finfo, correction, [('session_info', 'session_info')]),
                      (correction, contrasts, [('session_info', 'session_info')]),
                      (correction, l1d, [('session_info', 'session_info')]),
                      (contrasts, l1d, [('contrasts', 'contrasts')]),
                      (l1d, feat, [('fsf_files', 'fsf_file')]),
                      ])

