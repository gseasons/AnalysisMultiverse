#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:29:04 2021

@author: grahamseasons
"""
from updated.normalization.spatial_normalization import spatial_normalization
from nipype import Workflow, Node, IdentityInterface, Function
from nipype.interfaces.fsl import ExtractROI, MCFLIRT, SliceTimer, BET, FAST, FilterRegressor, GLM, UnaryMaths
from nipype.algorithms.rapidart import ArtifactDetect
import os
from updated.preprocessing.functions import get_wm, get_sink
from niworkflows.anat.ants import init_brain_extraction_wf
from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf



class preprocess(spatial_normalization):
    def __init__(self, task, pipeline, base_dir):
        super().__init__()
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        
    def construct(self, func_dic):
        preprocess = Workflow('preprocess')
        preprocess.base_dir = os.getcwd()
        inputnode = Node(IdentityInterface(fields=['bold', 'T1w', 'TR', 'mask']), name='inputnode')
        
        intermediates = ['smoothed', 'segmentations', 'warp', 'warp_field', 'brain', 'outliers', 'plots', 'invwarp', 'coregmat']
        
        self.preprocessing(preprocess, func_dic)
        self.get_warps(preprocess)
        
        preprocess.connect([(inputnode, preprocess.get_node('Fregistration'), [('T1w', 'T1w'),
                                                                              ('mask', 'mask')]),
                            (inputnode, preprocess.get_node('bet'), [('T1w', 'inputnode.in_files')]),
                            (inputnode, preprocess.get_node('slicetimer'), [('TR', 'time_repetition')]),
                            (inputnode, preprocess.get_node('extract'), [('bold', 'in_file')]),
                            (inputnode, preprocess.get_node('warp'), [('mask', 'ref_file')]),
                            (inputnode, preprocess.get_node('prelim'), [('mask', 'reference')]),
                            (inputnode, preprocess.get_node('dilateref'), [('mask', 'in_file')]),
                            (inputnode, preprocess.get_node('Fmni'), [('mask', 'mniMask')]),
                            (preprocess.get_node('bet_strip'), preprocess.get_node('warp'), [('out_file', 'in_file')]),
                            (preprocess.get_node('bet_strip'), preprocess.get_node('invwarp'), [('out_file', 'brain')]),
                            (preprocess.get_node('Fregistration'), preprocess.get_node('invwarp'), [('out_mat', 'coregmat')]),
                            (preprocess.get_node('Fmni'), preprocess.get_node('invwarp'), [('brainmask', 'brainmask')]),
                            (preprocess.get_node('bet_strip'), preprocess.get_node('prelim'), [('out_file', 'in_file')]),
                            (preprocess.get_node('bet_strip'), preprocess.get_node('decision'), [('out_file', 'mask')]),
                            (preprocess.get_node('bet_strip'), preprocess.get_node('dilatebrain'), [('out_file', 'in_file')]),
                            (preprocess.get_node('warp'), preprocess.get_node('Fmni'), [('field_file', 'warp')]),
                            ])
        
        outnode = Node(IdentityInterface(fields=['smoothed', 'segmentations', 'warp_file', 'outliers', 'brain', 'brainmask', 'unsmoothed', 'invwarp', 'keepreg', 'keepsmooth']), name='outnode')
        
        preprocess.connect([(preprocess.get_node('bet_strip'), outnode, [('out_file', 'brain')]),
                            (preprocess.get_node('Fmni'), outnode, [('segmentations', 'segmentations')]),
                            (preprocess.get_node('Fsmooth'), outnode, [('smooth', 'smoothed')]),
                            (preprocess.get_node('Fsmooth'), outnode, [('files', 'keepsmooth')]),
                            #(preprocess.get_node('warp'), outnode, [('field_file', 'warp_file')]),
                            (preprocess.get_node('Fmni'), outnode, [('warp', 'warp_file')]),
                            (preprocess.get_node('invwarp'), outnode, [('invwarp', 'invwarp')]),
                            #(preprocess.get_node('mcflirt'), outnode, [('par_file', 'mc_par')]),
                            (preprocess.get_node('art'), outnode, [('outlier_files', 'outliers')]),
                            #(preprocess.get_node('Fregistration'), outnode, [('out_mat', 'coregmat')]),
                            (preprocess.get_node('Fregistration'), outnode, [('files', 'keepreg')]),
                            (preprocess.get_node('fillmask'), outnode, [('out_file', 'brainmask')]),
                            ])
        
        if 'rest' in self.task:
            preprocess.connect([(preprocess.get_node('Fregress'), outnode, [('warped', 'unsmoothed')]),
                                ])
        else:
            preprocess.connect([(preprocess.get_node('Fmni'), outnode, [('warped', 'unsmoothed')]),
                                ])
        
        
        write = Node(Function(input_names=['base_dir', 'pipeline_st', 'task']+intermediates), name='write')
        write.inputs.function_str = get_sink(intermediates)
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        preprocess.connect([(preprocess.get_node('bet_strip'), write, [('out_file', 'brain')]),
                            (preprocess.get_node('Fsmooth'), write, [('smooth', 'smoothed')]),
                            (preprocess.get_node('Fregistration'), write, [('out_mat', 'coregmat')]),
                            (preprocess.get_node('fast'), write, [('tissue_class_files', 'segmentations')]),
                            (preprocess.get_node('warp'), write, [('field_file', 'warp_field'),
                                                                  ('warped_file', 'warp')]),
                            (preprocess.get_node('invwarp'), write, [('invwarp', 'invwarp')]),
                            (preprocess.get_node('art'), write, [('outlier_files', 'outliers'),
                                                                 ('plot_files', 'plots')]),
                            ])
        
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
        
        decision = Node(Function(input_names=['mask', 'mc_mean', 'mc', 'st', 'slice_correct', 'mean_vol'],
                                 output_names=['start_img', 'corrected_img', 'mask'], function=decision), name='decision')
        decision.inputs.mean_vol = ''
        decision.inputs.st = ''
        #NOTE: RIGHT NOW A WARP TO STANDARD SPACE MESSES THINGS UP WITH SEGMENTATION, ETC. -> PUT TOSTD IN OWN NODE, WARP SEGMENTATIONS, BRAINMASK, BOLD
        art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=3,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  plot_type='svg'),
                   name="art")
        
        fillmask = Node(UnaryMaths(operation='fillh'), name='fillmask')
        
        if 'rest' in self.task:
            func_str, input_names = function_str('regress', func_dic)
            Fregress = Node(Function(input_names=input_names,
                                     output_names=['warped']), name='Fregress')
            Fregress.inputs.function_str = func_str
            
            flow.connect([(flow.get_node('Fmni'), Fregress, [('warped', 'unsmoothed')]),
                          (flow.get_node('Fmni'), Fregress, [('brainmask', 'mask')]),
                          (flow.get_node('Fmni'), Fregress, [('segmentations', 'segmentations')]),
                          (mcflirt, Fregress, [('par_file', 'mc_par')]),
                          (Fregress, Fsmooth, [('warped', 'warped')]),
                          (fillmask, Fsmooth, [('out_file', 'mask')]),
                          ])
        else:
            flow.connect([(flow.get_node('Fmni'), Fsmooth, [('warped', 'warped')]),
                          (fillmask, Fsmooth, [('out_file', 'mask')]),
                          ])
        
        flow.connect([(extract, mcflirt, [('roi_file', 'in_file')]),
                      (mcflirt, slicetimer, [('out_file', 'in_file')]),
                      (mcflirt, decision, [('mean_img', 'mean_vol'),
                                           ('out_file', 'mc')]),
                      (mcflirt, art, [('par_file', 'realignment_parameters')]),
                      (slicetimer, decision, [('slice_time_corrected_file', 'st')]),
                      (decision, flow.get_node('Fregistration'), [('start_img', 'start_img'),
                                                                  ('corrected_img', 'corrected_img'),
                                                                  ('mask', 'brainmask')]),
                      (decision, flow.get_node('Fmni'), [('mask', 'brainmask')]),
                      (decision, flow.get_node('boldmask'), [('start_img', 'inputnode.in_file')]),
                      (flow.get_node('boldmask'), fillmask, [('outputnode.skull_stripped_file', 'in_file')]),
                      (flow.get_node('Fregistration'), art, [('warped', 'realigned_files')]),
                      ])
        
    def coregistration(self, flow, func_dic):
        from updated.preprocessing.functions import function_str, strip_container
        bet = init_brain_extraction_wf(name='bet')
        bet_strip = Node(Function(input_names='in_file', output_names='out_file', function=strip_container), name='bet_strip')
        #bet = Node(BET(output_type='NIFTI_GZ'), name='bet')
        fast = Node(FAST(output_type='NIFTI_GZ', segments=True), name='fast')
        
        func_str, input_names = function_str('registration', func_dic)
        
        Fregistration = Node(Function(input_names=input_names,
                                      output_names=['out_mat', 'warped', 'files']), 
                             name='Fregistration')
        Fregistration.inputs.function_str = func_str
        
        func_str, input_names = function_str('mni', func_dic)
        
        Fmni = Node(Function(input_names=input_names,
                             output_names=['warped', 'brainmask', 'segmentations', 'warp']), 
                             name='Fmni')
        Fmni.inputs.function_str = func_str
        
        boldmask = init_enhance_and_skullstrip_bold_wf(name='boldmask', pre_mask=True)
        
        flow.connect([(bet, bet_strip, [('outputnode.out_file', 'in_file')]),
                      (bet_strip, fast, [('out_file', 'in_files')]),
                      (bet_strip, Fregistration, [('out_file', 'bet')]),
                      (bet_strip, Fmni, [('out_file', 'brain')]),
                      (fast, Fregistration, [(('tissue_class_files', get_wm), 'wm_file')]),
                      (fast, Fmni, [('tissue_class_files', 'segmentations')]),
                      (Fregistration, Fmni, [('warped', 'warped')]),
                      (Fregistration, Fmni, [('out_mat', 'out_mat')]),
                      (Fmni, boldmask, [('brainmask', 'inputnode.pre_mask')]),
                      ])
