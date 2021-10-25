#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 12:26:40 2021

@author: grahamseasons
"""
from nipype import Workflow, Node, MapNode, Function, IdentityInterface
from nipype.interfaces.fsl import L2Model, FLAMEO, Merge
import os

def groupscans(copes, varcopes):
    outcopes = []
    outvarcopes = []
    multiple = type(copes[0])
    if multiple == list:
        numcon = len(copes[0])
    else:
        return [copes], [varcopes], len(copes)
        
    for i in range(numcon):
        coperuns = [run[i] for run in copes]
        varcoperuns = [run[i] for run in varcopes]
        outcopes.append(coperuns)
        outvarcopes.append(varcoperuns)
        
    return outcopes, outvarcopes, numcon
            

class level2:
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
        
        
    def l2(self, flow):
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
        
        
        
        