#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:00:22 2021

@author: grahamseasons
"""
from nipype import Workflow, Node, MapNode, Function, IdentityInterface
from nipype.interfaces.fsl import MultipleRegressDesign, FLAMEO, Merge
import os
    
class level3:
    def __init__(self, task, pipeline, base_dir):
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self):
        from l3_analysis.functions import get_sink, remove_container, mniMask
        level3 = Workflow('level3')
        level3.base_dir = os.getcwd()
        #remove container from outputs
        inputnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'mask', 'covariates', 'subjects']), name='inputnode')
        
        intermediates = ['copes', 'varcopes', 'flameo_stats']
        files = []#stats/file
        
        self.l3(level3)
        
        level3.connect([(inputnode, level3.get_node('group'), [('covariates', 'covariate'),
                                                               ('subjects', 'subjects')]),
                        (inputnode, level3.get_node('groupcontrast'), [('copes', 'copes'),
                                                                       ('varcopes', 'varcopes')]),
                        (inputnode, level3.get_node('flameo'), [(('mask', mniMask), 'mask_file')]),
                        ])
        
        outnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'zstats', 'flameo_stats', 'mask']), name='outnode')
        
        level3.connect([(level3.get_node('flameo'), outnode, [(('copes', remove_container), 'copes'),
                                                              (('var_copes', remove_container), 'varcopes'),
                                                              (('zstats', remove_container), 'zstats'),
                                                              (('stats_dir', remove_container), 'flameo_stats')]),
                        (inputnode, outnode, [(('mask', mniMask), 'mask')]),
                        ])
        
        write = Node(Function(input_names=['base_dir', 'pipeline_st', 'task']+intermediates), name='write')
        write.inputs.function_str = get_sink(intermediates, files)
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        level3.connect([(level3.get_node('flameo'), write, [(('copes', remove_container), 'copes'),
                                                            (('var_copes', remove_container), 'varcopes'),
                                                            (('zstats', remove_container), 'zstats'),
                                                            (('stats_dir', remove_container), 'flameo_stats')]),
                        ])
        
        return level3
        
        
    def l3(self, flow):
        from l3_analysis.functions import group_contrast
        self.setup(flow)
        
        l3model = Node(MultipleRegressDesign(), name='l3model')
        
        groupcontrast = Node(Function(input_names=['copes', 'varcopes'],
                                  output_names=['copes', 'varcopes'], 
                                  function=group_contrast), name='groupcontrast')
        
        mergecopes = MapNode(Merge(dimension='t'), iterfield='in_files', name='mergecopes')
        mergevarcopes = MapNode(Merge(dimension='t'), iterfield='in_files', name='mergevarcopes')
        
        flameo = MapNode(FLAMEO(), iterfield=['cope_file', 'var_cope_file'], name='flameo')
        
        flow.connect([(flow.get_node('group'), l3model, [('EVs', 'regressors'),
                                                         ('contrasts', 'contrasts'),
                                                         ('group_ids', 'groups')]),
                      (groupcontrast, mergecopes, [('copes', 'in_files')]),
                      (groupcontrast, mergevarcopes, [('varcopes', 'in_files')]),
                      (mergecopes, flameo, [('merged_file', 'cope_file')]),
                      (mergevarcopes, flameo, [('merged_file', 'var_cope_file')]),
                      (l3model, flameo, [('design_con', 't_con_file'),
                                         ('design_grp', 'cov_split_file'),
                                         ('design_fts', 'f_con_file'),
                                         ('design_mat', 'design_file')]),
                      ])
        
        
    def setup(self, flow):
        from l3_analysis.functions import t_test
        
        group = Node(Function(input_names=['covariate', 'subjects'],
                              output_names=['EVs', 'contrasts', 'group_ids'], function=t_test), name='group')
        
        flow.add_nodes([group])
        
        
        
        
        
        
        
        
        
        