#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:31:44 2021

@author: grahamseasons
"""
from nipype import Node, Workflow, IdentityInterface, Function
from updated.split_half.functions import get_sink
import os

class split:
    def __init__(self, task, pipeline, base_dir):
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self):
        split = Workflow('split')
        split.base_dir = os.getcwd()
        
        inputnode = Node(IdentityInterface(fields=['corrected', 'mask', 'bold', 'covariates']), name='inputnode')
        
        self.prediction(split)
        self.reproducibility(split)
        self.distance(split)
        
        split.connect([(inputnode, split.get_node('repro'), [('corrected', 'stats'),
                                                             ('mask', 'mask')]),
                       (inputnode, split.get_node('predict'), [('bold', 'bold'),
                                                               ('covariates', 'covariate'),
                                                               ('mask', 'mask')]),
                       ])
        
        intermediates = ['P_txt', 'R_txt', 'score_txt']
        
        write = Node(Function(input_names=['base_dir', 'pipeline_st', 'task']+intermediates), name='write')
        write.inputs.function_str = get_sink(intermediates)
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        split.connect([(split.get_node('predict'), write, [('P_txt', 'P_txt')]),
                       (split.get_node('repro'), write, [('R_txt', 'R_txt')]),
                       (split.get_node('calculation'), write, [('score_txt', 'score_txt')]),
                       ])
        
        return split
        
    def distance(self, flow):
        from updated.split_half.functions import calc
        calculation = Node(Function(input_names=['R', 'P'],
                                    output_names=['score', 'score_txt'], function=calc), name='calculation')
        
        flow.connect([(flow.get_node('repro'), calculation, [('R', 'R')]),
                      (flow.get_node('predict'), calculation, [('P', 'P')]),
                      ])
        
    def reproducibility(self, flow):#INPUT WILL BE LEVEL3 outputs, JOINED POST GROUP SPLITTING i.e. a bunch of pairs
        from updated.split_half.functions import compare
        
        repro = Node(Function(input_names=['stats', 'mask'],
                              output_names=['R', 'R_txt'], function=compare), name='repro')
        
        flow.add_nodes([repro])
        
    def prediction(self, flow):
        from updated.split_half.functions import predict
        
        predict = Node(Function(input_names=['covariate', 'bold', 'mask'],
                                output_names=['P', 'P_txt'], function=predict), name='predict')
        
        flow.add_nodes([predict])