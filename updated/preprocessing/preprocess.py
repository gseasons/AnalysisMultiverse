#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:29:04 2021

@author: grahamseasons
"""
from updated.normalization.spatial_normalization import spatial_normalization
from nipype import Workflow, Node, IdentityInterface, JoinNode, MapNode, Function
from nipype.interfaces.fsl import ExtractROI, MCFLIRT, SliceTimer, BET, FAST
from nipype.algorithms.rapidart import ArtifactDetect
import os
from updated.functions import traverse
from updated.preprocessing.functions import get_wm

class preprocess(spatial_normalization):
    def __init__(self, task, pipeline, base_dir):
        super().__init__()
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self, param_dic, func_dic):
        from updated.preprocessing.functions import out
        preprocess = Workflow('preprocess')
        preprocess.base_dir = os.getcwd()
        suffix = '_pre'
        inputnode = Node(IdentityInterface(fields=['bold', 'T1w', 'TR', 'mask']), name='inputnode')
        out_fields = ['smoothed', 'segmentations', 'warp_file', 'outliers', 'mc_par', 'brain']
        
        self.preprocessing(preprocess, func_dic)
        self.get_warps(preprocess)
        
        preprocess.connect([(inputnode, preprocess.get_node('Fregistration'), [('T1w', 'T1w'),
                                                                              ('mask', 'mask')]),
                            (inputnode, preprocess.get_node('bet'), [('T1w', 'in_file')]),
                            (inputnode, preprocess.get_node('slicetimer'), [('TR', 'time_repetition')]),
                            (inputnode, preprocess.get_node('extract'), [('bold', 'in_file')]),
                            (inputnode, preprocess.get_node('warp'), [('mask', 'ref_file')]),
                            (inputnode, preprocess.get_node('prelim'), [('mask', 'reference')]),
                            (inputnode, preprocess.get_node('dilateref'), [('mask', 'in_file')]),
                            (preprocess.get_node('bet'), preprocess.get_node('warp'), [('out_file', 'in_file')]),
                            (preprocess.get_node('bet'), preprocess.get_node('prelim'), [('out_file', 'in_file')]),
                            (preprocess.get_node('bet'), preprocess.get_node('dilatebrain'), [('out_file', 'in_file')]),
                            (preprocess.get_node('warp'), preprocess.get_node('Fregistration'), [('field_file', 'warp_file')]),
                            ])
        
        
        
        joinsource = traverse(param_dic, preprocess, suffix)
        
        outnode = Node(IdentityInterface(fields=out_fields), name='outnode')
        
        if joinsource:
            outnode = JoinNode(IdentityInterface(fields=['smoothed', 'segmentations', 'warp_file', 'outliers', 'mc_par', 'brain']),
                               name='outnode', joinsource=joinsource, joinfield=['smoothed', 'segmentations', 'warp_file', 'outliers', 'mc_par', 'brain'])
        else:
            outnode = Node(IdentityInterface(fields=['smoothed', 'segmentations', 'warp_file', 'outliers', 'mc_par', 'brain']), name='outnode')
        
        preprocess.connect([(preprocess.get_node('bet'), outnode, [('out_file', 'brain')]),
                            (preprocess.get_node('fast'), outnode, [('partial_volume_files', 'segmentations')]),
                            (preprocess.get_node('Fsmooth'), outnode, [('smooth', 'smoothed')]),
                            (preprocess.get_node('warp'), outnode, [('field_file', 'warp_file')]),
                            (preprocess.get_node('mcflirt'), outnode, [('par_file', 'mc_par')]),
                            (preprocess.get_node('art'), outnode, [('outlier_files', 'outliers')]),
                            ])
        
# =============================================================================
#         write = MapNode(Function(input_names=['base_dir', 'pipeline_st', 'task', 'param_dic', 'slicetiming', 'smooth', 'mc_par', 'mc_img', 'mc_mean', 'brain', 'outliers', 'plots', 'out_mat', 'warped', 'segment', 'warp_file'],
#                                  function=out), iterfield=['slicetiming', 'smooth', 'mc_par', 'mc_img', 'mc_mean', 'brain', 'outliers', 'plots', 'out_mat', 'warped', 'segment', 'warp_file'],
#                         name='write')
#         
#         write.inputs.base_dir = self.base_dir
#         write.inputs.pipeline_st = self.pipeline
#         write.inputs.task = self.task
# =============================================================================
        
        
        
        return preprocess
    
    def preprocessing(self, flow, func_dic):
        from updated.preprocessing.functions import function_str, decision
        self.coregistration(flow, func_dic)
        
        extract = Node(ExtractROI(t_size=-1, output_type='NIFTI_GZ'), name='extract')
        mcflirt = Node(MCFLIRT(save_plots=True, output_type='NIFTI_GZ'), name='mcflirt')
        #TURN ON/OFF
        slicetimer = Node(SliceTimer(index_dir=False, interleaved=True, output_type='NIFTI_GZ'), name='slicetimer')
        
        func_str, input_names = function_str('smooth', func_dic)
        Fsmooth = Node(Function(input_names=input_names,
                                   output_names=['smooth', 'files']), name='Fsmooth')
        Fsmooth.inputs.function_str = func_str
        
        decision = Node(Function(input_names=['mc_mean', 'mc', 'st', 'slice_correct', 'mean_vol'],
                                 output_names=['start_img', 'corrected_img'], function=decision), name='decision')
        decision.inputs.mean_vol = ''
        decision.inputs.st = ''
        
        art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=3,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  plot_type='svg'),
                   name="art")
        
        flow.connect([(extract, mcflirt, [('roi_file', 'in_file')]),
                      (mcflirt, slicetimer, [('out_file', 'in_file')]),
                      (mcflirt, decision, [('mean_img', 'mean_vol'),
                                           ('out_file', 'mc')]),
                      (mcflirt, art, [('par_file', 'realignment_parameters')]),
                      (slicetimer, decision, [('slice_time_corrected_file', 'st')]),
                      (decision, flow.get_node('Fregistration'), [('start_img', 'start_img'),
                                                                 ('corrected_img', 'corrected_img')]),
                      (flow.get_node('Fregistration'), Fsmooth, [('warped', 'warped')]),
                      (flow.get_node('Fregistration'), art, [('warped', 'realigned_files')]),
                      ])
        
    def coregistration(self, flow, func_dic):
        from updated.preprocessing.functions import function_str
        bet = Node(BET(output_type='NIFTI_GZ'), name='bet')
        fast = Node(FAST(output_type='NIFTI_GZ'), name='fast')
        
        func_str, input_names = function_str('registration', func_dic)
        
        Fregistration = Node(Function(input_names=input_names,
                                output_names=['out_mat', 'warped', 'files']), 
                                name='Fregistration')
        Fregistration.inputs.function_str = func_str
        
        flow.connect([(bet, fast, [('out_file', 'in_files')]),
                      (bet, Fregistration, [('out_file', 'bet')]),
                      (fast, Fregistration, [(('partial_volume_files', get_wm), 'wm_file')]),
                      ])
