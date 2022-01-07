#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 12:35:09 2021

@author: grahamseasons
"""
from nipype import Node, Workflow, Function, IdentityInterface
import os

class correction:
    def __init__(self, task, pipeline, base_dir):
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self):
        from correction.functions import get_sink
        correction = Workflow('correction')
        correction.base_dir = os.getcwd()
        
        inputnode = Node(IdentityInterface(fields=['zstat', 'copes', 'mask']), name='inputnode')
        
        intermediates = ['corrected_files']
        
        self.decision(correction)
        
        correction.connect([(inputnode, correction.get_node('correct'), [('zstat', 'zstat'),
                                                                         ('copes', 'copes'),
                                                                         ('mask', 'mask')]),
                            ])
        
        outnode = Node(IdentityInterface(fields=['corrected_files']), name='outnode')
        
        correction.connect([(correction.get_node('correct'), outnode, [('corrected_files', 'corrected_files')])])
        
        write = Node(Function(input_names=['base_dir', 'pipeline_st', 'task']+intermediates), name='write')
        write.inputs.function_str = get_sink(intermediates)
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        correction.connect([(outnode, write, [('corrected_files', 'corrected_files')])])
        
        return correction
        
        
    def decision(self, flow):
        from correction.functions import correction
        correct = Node(Function(input_names=['zstat', 'copes', 'mask', 'cor'],
                                output_names=['corrected_files'], 
                                function=correction), name='correct')
        
        flow.add_nodes([correct])