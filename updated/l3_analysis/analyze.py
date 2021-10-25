#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:00:22 2021

@author: grahamseasons
"""
from nipype import Workflow, Node, MapNode, Function, IdentityInterface
from nipype.interfaces.fsl import MultipleRegressDesign, FLAMEO, Merge
import os

def group_contrast(lst):
    subs = len(lst)
    contrasts = len(lst[0])
    grouped = []
    for con in range(contrasts):
        grouped.append([lst[sub][con] for sub in range(subs)])
        
    return grouped
    
class level3:
    def __init__(self, task, pipeline, base_dir):
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self, func_dic):
        level2 = Workflow('level1')
        level2.base_dir = os.getcwd()
        
        inputnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'mask']), name='inputnode')
        
        self.l2(level2)
        
        level2.connect([(inputnode, level2.get_node('groupruns'), [('copes', 'copes'),
                                                                   ('varcopes', 'varcopes')]),
                        (inputnode, level2.get_node('flameo'), [('mask', 'mask')]),
                        ])
        
        outnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'flameo_stats']))
        
        level2.connect([(level2.get_node('flameo'), outnode, [('copes', 'copes'),
                                                              ('varcopes', 'varcopes')]),
                        (level2.get_node('flameo'), outnode, [('flameo_stats', 'flameo_stats')]),
                        ])
        
        
    def l3(self, flow):
        l3model = Node(MultipleRegressDesign(), name='l3model')
        
        groupruns = Node(Function(input_names=['copes', 'varcopes'],
                                  output_names=['copes', 'varcopes', 'num_copes'], 
                                  function=group_contrast), name='groupruns')
        
        mergecopes = MapNode(Merge(dimension='t'), iterfield='in_files', name='mergecopes')
        mergevarcopes = MapNode(Merge(dimension='t'), iterfield='in_files', name='mergevarcopes')
        
        flameo = MapNode(FLAMEO(), iterfield=['cope_file', 'var_cope_file'], name='flameo')
        
        
        flow.connect([(flow.get_node('group'), l3model, [('EVs', 'regressors'),
                                                         ('contrasts', 'contrasts'),
                                                         ('group_ids', 'groups')]),
                      ])
        
        flow.connect([(groupruns, mergecopes, [('copes', 'in_files')]),
                      (groupruns, mergevarcopes, [('varcopes', 'in_files')]),
                      (groupruns, l2model, [('num_copes', 'num_copes')]),
                      (l2model, flameo, [('design_mat', 'design_mat'),
                                         ('design_con', 'design_con'),
                                         ('design_grp', 'design_grp')]),
                      (mergecopes, flameo, [('merged_file', 'cope_file')]),
                      (mergevarcopes, flameo, [('merged_file', 'var_cope_file')]),
                      ])
        
        
    def setup(self, flow):
        from updated.l3_analysis.functions import t_test
        
        group = Node(Function(input_names=['covariates', 'subjects'],
                              output_names=['EVs', 'contrasts', 'group_ids'], function=t_test), name='group')
        
        flow.add_nodes([group])
        
        
        
        
        
        
        
        
        
        