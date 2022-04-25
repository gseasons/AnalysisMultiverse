#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:29:04 2021

@author: grahamseasons
"""
from normalization.spatial_normalization import spatial_normalization
from nipype import Workflow, Node, IdentityInterface, Function, DataSink
from nipype.interfaces.fsl import ExtractROI, MCFLIRT, SliceTimer, FAST, UnaryMaths
from nipype.algorithms.rapidart import ArtifactDetect
import os, glob
from preprocessing.functions import get_wm, get_sink
#from niworkflows.anat.ants import init_brain_extraction_wf
from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf

class preprocess(spatial_normalization):
    def __init__(self, task, pipeline, base_dir, data_dir):
        super().__init__()
        self.task = task
        self.pipeline = pipeline
        self.base_dir = base_dir
        self.data_dir = data_dir
        
    def construct(self, func_dic, subjects):
        preprocess = Workflow('preprocess')
        preprocess.base_dir = os.getcwd()
        inputnode = Node(IdentityInterface(fields=['bold', 'T1w', 'TR', 'mask']), name='inputnode')
        
        intermediates = ['smoothed', 'segmentations', 'warp', 'warp_field', 'brain', 'outliers', 'plots', 'invwarp', 'coregmat']
        
        brain_extracted = self.preprocessing(preprocess, func_dic, subjects)
        self.get_warps(preprocess)
        
        preprocess.connect([(inputnode, preprocess.get_node('Fregistration'), [('T1w', 'T1w')]),
                                                                              #('mask', 'mask')]),
                            #(inputnode, preprocess.get_node('bet'), [('T1w', 'inputnode.in_files')]),
                            #(inputnode, preprocess.get_node('bet'), [('T1w', 'T1w')]),
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
        
        if brain_extracted:
            preprocess.connect(inputnode, 'T1w', preprocess.get_node('bet'), 'T1w')
        else:
            preprocess.connect(inputnode, 'T1w', preprocess.get_node('bet'), 'inputnode.in_files')
        
        outnode = Node(IdentityInterface(fields=['smoothed', 'segmentations', 'warp_file', 'brain', 'brainmask', 'outliers', 'unsmoothed', 'invwarp']), name='outnode')
        
        preprocess.connect([(preprocess.get_node('bet_strip'), outnode, [('out_file', 'brain')]),
                            (preprocess.get_node('Fmni'), outnode, [('segmentations', 'segmentations')]),
                            (preprocess.get_node('Fsmooth'), outnode, [('smooth', 'smoothed')]),
                            #(preprocess.get_node('Fsmooth'), outnode, [('files', 'keepsmooth')]),
                            #(preprocess.get_node('warp'), outnode, [('field_file', 'warp_file')]),
                            (preprocess.get_node('Fmni'), outnode, [('warp', 'warp_file')]),
                            (preprocess.get_node('invwarp'), outnode, [('invwarp', 'invwarp')]),
                            #(preprocess.get_node('mcflirt'), outnode, [('par_file', 'mc_par')]),
                            (preprocess.get_node('art'), outnode, [('outlier_files', 'outliers')]),
                            #(preprocess.get_node('Fregistration'), outnode, [('out_mat', 'coregmat')]),
                            #(preprocess.get_node('Fregistration'), outnode, [('files', 'keepreg')]),
                            (preprocess.get_node('fillmask'), outnode, [('out_file', 'brainmask')]),
                            ])
        
        if 'rest' in self.task:
            preprocess.connect([(preprocess.get_node('Fregress'), outnode, [('forreho', 'unsmoothed')]),
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
    
    def preprocessing(self, flow, func_dic, subjects):
        from preprocessing.functions import function_str, decision
        brain_extracted = self.coregistration(flow, func_dic, subjects)
        
        extract = Node(ExtractROI(t_size=-1, output_type='NIFTI_GZ'), name='extract', n_procs=3, mem_gb=0.5)
        mcflirt = Node(MCFLIRT(save_plots=True, output_type='NIFTI_GZ'), name='mcflirt', n_procs=2, mem_gb=0.7)
        
        slicetimer = Node(SliceTimer(index_dir=False, interleaved=True, output_type='NIFTI_GZ'), name='slicetimer', n_procs=2, mem_gb=0.4)
            
        func_str, input_names = function_str('smooth', func_dic)
        Fsmooth = Node(Function(input_names=input_names,
                                output_names=['smooth']), name='Fsmooth', n_procs=3, mem_gb=2)
        Fsmooth.inputs.function_str = func_str
        
        decision = Node(Function(input_names=['mask', 'mc_mean', 'mc', 'st', 'slice_correct', 'mean_vol'],
                                 output_names=['start_img', 'corrected_img', 'mask'], function=decision), name='decision', mem_gb=0.3)
        decision.inputs.mean_vol = ''
        decision.inputs.st = ''
        
        art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=3,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  plot_type='svg'),
                   name="art", mem_gb=0.5)
        
        fillmask = Node(UnaryMaths(operation='fillh'), name='fillmask')
        
        func_str, input_names = function_str('regress', func_dic)
        Fregress = Node(Function(input_names=input_names,
                                 output_names=['warped', 'forreho']), name='Fregress', n_procs=5, mem_gb=6)
        Fregress.inputs.function_str = func_str
        
        if 'rest' in self.task:
            Fregress.inputs.rest = True
        else:
            Fregress.inputs.rest = False
            
        flow.connect([(extract, mcflirt, [('roi_file', 'in_file')]),
                      (mcflirt, slicetimer, [('out_file', 'in_file')]),
                      (mcflirt, decision, [('mean_img', 'mean_vol'),
                                           ('out_file', 'mc')]),
                      (mcflirt, art, [('par_file', 'realignment_parameters')]),
                      (slicetimer, decision, [('slice_time_corrected_file', 'st')]),
                      (decision, flow.get_node('Fregistration'), [('start_img', 'start_img'),
                                                                  ('corrected_img', 'corrected_img')]),
                             #                                     ('mask', 'brainmask')]), SHOULDN'T BE NEEDED
                      (decision, flow.get_node('Fmni'), [('mask', 'brainmask')]),
                      (decision, flow.get_node('Fmni'), [('start_img', 'start_img')]),
                      (flow.get_node('Fmni'), flow.get_node('boldmask'), [('start_img', 'inputnode.in_file')]),
                      #(decision, flow.get_node('boldmask'), [('start_img', 'inputnode.in_file')]),
                      (flow.get_node('boldmask'), fillmask, [('outputnode.skull_stripped_file', 'in_file')]),
                      (flow.get_node('Fregistration'), art, [('warped', 'realigned_files')]),
                      (flow.get_node('Fmni'), Fregress, [('warped', 'unsmoothed')]),
                      #MAYBE REPLACE THIS WITH OUTPUT OF FILLMASK
                      #(flow.get_node('fillmask'), Fregress, [('out_file', 'mask')]),
                      (flow.get_node('Fmni'), Fregress, [('brainmask', 'mask')]),
                      (flow.get_node('Fmni'), Fregress, [('segmentations', 'segmentations')]),
                      #(art, Fregress, [('outlier_files', 'outliers')]),
                      (mcflirt, Fregress, [('par_file', 'mc_par')]),
                      (Fregress, Fsmooth, [('warped', 'warped')]),
                      (fillmask, Fsmooth, [('out_file', 'mask')]),
                      ])
        
        return brain_extracted
        
    def coregistration(self, flow, func_dic, subjects):
        from preprocessing.functions import function_str, strip_container
        from preprocessing.workflows import check4brains, betwrite
        
        bet = Node(Function(input_names=['data_dir', 'T1w'], output_names='out_file', function=check4brains), name='bet')
        bet.inputs.data_dir = self.data_dir
        brain_extracted = True
        
        for sub in subjects:
            masked = glob.glob(self.data_dir + '/brain_extracted/_subject_{sID}/*masked.nii.gz'.format(sID=sub))
            if not masked:
                from niworkflows.anat.ants import init_brain_extraction_wf
                bet = init_brain_extraction_wf(name='bet', omp_nthreads=2, mem_gb=3)
                brain_extracted = False
                break
        
        bet_strip = Node(Function(input_names='in_file', output_names='out_file', function=strip_container), name='bet_strip')
        fast = Node(FAST(output_type='NIFTI_GZ', segments=True), name='fast', n_procs=2, mem_gb=1.7)
        
        if sub != subjects[-1]:
            betwrite = Node(Function(input_names=['data_dir', 'bet_'], function=betwrite), name='betwrite')
            betwrite.inputs.data_dir = self.data_dir
            flow.connect(bet, 'outputnode.out_file', bet_strip, 'in_file')
            flow.connect(bet_strip, 'out_file', betwrite, 'bet_')
        else:
            flow.connect(bet, 'out_file', bet_strip, 'in_file')
        
        func_str, input_names = function_str('registration', func_dic)
        
        Fregistration = Node(Function(input_names=input_names,
                                      output_names=['out_mat', 'warped']), 
                             name='Fregistration', n_procs=3, mem_gb=0.9)
        Fregistration.inputs.function_str = func_str
        
        func_str, input_names = function_str('mni', func_dic)
        
        Fmni = Node(Function(input_names=input_names,
                             output_names=['warped', 'brainmask', 'segmentations', 'warp', 'start_img']), 
                             name='Fmni', n_procs=4, mem_gb=1.5)
        Fmni.inputs.function_str = func_str
        
        boldmask = init_enhance_and_skullstrip_bold_wf(name='boldmask', pre_mask=True)
        
        flow.connect([(bet_strip, fast, [('out_file', 'in_files')]),
                      (bet_strip, Fregistration, [('out_file', 'bet')]),
                      (bet_strip, Fmni, [('out_file', 'brain')]),
                      (fast, Fregistration, [(('tissue_class_files', get_wm), 'wm_file')]),
                      (fast, Fmni, [('tissue_class_files', 'segmentations')]),
                      (Fregistration, Fmni, [('warped', 'warped')]),
                      (Fregistration, Fmni, [('out_mat', 'out_mat')]),
                      (Fmni, boldmask, [('brainmask', 'inputnode.pre_mask')]),
                      ])
        
        return brain_extracted
