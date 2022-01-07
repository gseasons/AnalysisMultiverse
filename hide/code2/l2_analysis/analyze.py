#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 12:26:40 2021

@author: grahamseasons
"""
from nipype import Workflow, Node, MapNode, Function, IdentityInterface
from nipype.interfaces.fsl import L2Model, FLAMEO, Merge
import os            

class level2:
    def __init__(self, task, pipeline, base_dir):
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self):
        from l2_analysis.functions import get_sink
        level2 = Workflow('level2')
        level2.base_dir = os.getcwd()
        
        inputnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'mask']), name='inputnode')
        
        intermediates = ['copes', 'varcopes', 'flameo_stats']
        files = []#stats/file
        
        self.l2(level2)
        
        level2.connect([(inputnode, level2.get_node('groupruns'), [('copes', 'copes'),
                                                                   ('varcopes', 'varcopes')]),
                        (inputnode, level2.get_node('flameo'), [('mask', 'mask_file')]),
                        ])
        
        outnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'flameo_stats']), name='outnode')
        
        level2.connect([(level2.get_node('flameo'), outnode, [('copes', 'copes'),
                                                              ('var_copes', 'varcopes'),
                                                              ('zstats', 'zstats'),
                                                              ('stats_dir', 'flameo_stats')]),
                        ])
        
        write = Node(Function(input_names=['base_dir', 'pipeline_st', 'task']+intermediates), name='write')
        write.inputs.function_str = get_sink(intermediates, files)
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        level2.connect([(level2.get_node('flameo'), write, [('copes', 'copes'),
                                                            ('var_copes', 'varcopes'),
                                                            ('stats_dir', 'flameo_stats')]),
                        ])
        
        return level2
        
        
    def l2(self, flow):
        from l2_analysis.functions import groupscans
        groupruns = Node(Function(input_names=['copes', 'varcopes'],
                                  output_names=['copes', 'varcopes', 'num_copes'], 
                                  function=groupscans), name='groupruns')
        
        mergecopes = MapNode(Merge(dimension='t'), iterfield='in_files', name='mergecopes')
        mergevarcopes = MapNode(Merge(dimension='t'), iterfield='in_files', name='mergevarcopes')
        l2model = Node(L2Model(), name='l2model')
        flameo = MapNode(FLAMEO(), iterfield=['cope_file', 'var_cope_file'], name='flameo')
        
        flow.connect([(groupruns, mergecopes, [('copes', 'in_files')]),
                      (groupruns, mergevarcopes, [('varcopes', 'in_files')]),
                      (groupruns, l2model, [('num_copes', 'num_copes')]),
                      (l2model, flameo, [('design_mat', 'design_mat'),
                                         ('design_con', 'design_con'),
                                         ('design_grp', 'design_grp')]),
                      (mergecopes, flameo, [('merged_file', 'cope_file')]),
                      (mergevarcopes, flameo, [('merged_file', 'var_cope_file')]),
                      ])
        
        
        
        