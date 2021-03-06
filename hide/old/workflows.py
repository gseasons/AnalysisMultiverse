#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 12:25:19 2021

@author: grahamseasons
"""
from os.path import join as opj
from bids.layout import BIDSLayout
from nipype.interfaces.io import BIDSDataGrabber
from nipype import Workflow, Node, JoinNode, Function, DataSink, IdentityInterface
from os.path import join as opj
import os
from nipype.interfaces.fsl import (BET, FAST, FLIRT, Threshold)
from nipype import Workflow, Node, MapNode

from os.path import join as opj
import os
import json
from nipype.interfaces.fsl import (BET, ExtractROI, FAST, FLIRT, ImageMaths, ImageStats,
                                   MCFLIRT, SliceTimer, Threshold, SUSAN, IsotropicSmooth)
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.algorithms.rapidart import ArtifactDetect
from nipype import Workflow, Node

from os.path import join as opj
from nipype.interfaces.fsl import (FEATModel, FEATRegister, FEAT, FILMGLS, GLM, Level1Design)
from nipype.algorithms.modelgen import SpecifyModel
from nipype import Node, Workflow, Function

from nipype.interfaces.fsl import FNIRT, ApplyWarp, FLIRT
from nipype import Node, MapNode, Workflow, SelectFiles, IdentityInterface
import os
from os.path import join as opj

from nipype import Node, JoinNode, MapNode, Workflow, SelectFiles, IdentityInterface, Function
from nipype.interfaces.fsl import L2Model, FLAMEO, Cluster, ImageMaths, Merge
import os
from os.path import join as opj
from functions import spatial_normalization

#mask_file = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain_mask.nii.gz')

#PREPROCESSING
class preprocess(spatial_normalization):
    def __init__(self, base_dir, task):
        super().__init__()
        self.base_dir = base_dir
        self.task = task
        self.pipeline = 0
        
    def construct(self, pre):
        preprocess = Workflow('preprocess')
        preprocess.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['bold', 'T1w', 'TR', 'mask']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['smoothed', 'outliers', 'plots', 'mc_par',
                                                 'warped_mean', 'warped', 'brain', 'reg_out_mat',
                                                 'warp_file', 'segmentations']), name='outnode')
        
        im = self.improcess_flow(pre)
        #reg = self.coreg_flow()
        
        preprocess.connect([(inputnode, im, [('bold', 'inputnode_im.bold'),
                                             ('T1w', 'inputnode_im.T1w'),
                                             ('TR', 'inputnode_im.TR'),
                                             #('frac_mask', 'inputnode_im.frac_mask'),
                                             #('discard', 'inputnode_im.discard'),
                                             #('dof_mc', 'inputnode_im.dof_mc'),
                                             #('fwhm', 'inputnode_im.fwhm'),
                                             #('cost_mc', 'inputnode_im.cost_mc'),
                                             #('susan', 'inputnode_im.susan'),
                                             #('warplater', 'inputnode_im.warplater'),
                                             #('warp_file', 'inputnode.warp_file'),
                                             ('mask', 'inputnode_im.mask')]),
# =============================================================================
#                             (inputnode, im, [('bet_frac', 'inputnode_im.bet_frac'),
#                                               ('robust', 'inputnode_im.robust'),
#                                               ('wm_thresh', 'inputnode_im.wm_thresh'),
#                                               ('dof_f', 'inputnode_im.dof_f'),
#                                               ('bbr_type', 'inputnode_im.bbr_type'),
#                                               ('interp', 'inputnode_im.interp'),
#                                               ('cost', 'inputnode_im.cost'),
#                                               ('bins', 'inputnode_im.bins'),
#                                               ('iso', 'inputnode_im.iso'),
#                                               ('bbr', 'inputnode_im.bbr')]),
# =============================================================================
                            (im, outnode, [('outnode.smoothed', 'smoothed'),
                                           ('outnode.outliers', 'outliers'),
                                           ('outnode.plots', 'plots'),
                                           ('outnode.mc_par', 'mc_par')]),
                            (im, outnode, [('outnode.warped_mean', 'warped_mean'),
                                            ('outnode.warped', 'warped'),
                                            ('outnode.brain', 'brain'),
                                            ('outnode.out_mat', 'reg_out_mat'),
                                            ('outnode.segmentations', 'segmentations'),
                                            ('outnode.warp_file', 'warp_file')]),
                                            #('outnode.invwarp', 'invwarp')]),
                            ])
        
        return preprocess
        
        
        #frac_mask=0.3, TR, options, susan=False, bbr=True, discard=4, dof_mc=6, fwhm=4, cost_mc='normcorr'
    def improcess_flow(self, pre):
        #ASSUMES INTERLEAVED SLICE TIMING, USES NIPYPE ERROR DETECTION INSTEAD OF MELODIC
        imageproc = Workflow(name='imageproc')
        imageproc.base_dir = os.getcwd() #self.base_dir
        im_var = pre['improcess']
        co_var = pre['coreg']
        
        rng = len(im_var['PIPELINE'][0][1])
        
        fields=['bold', 'T1w', 'susan', 'TR', 'frac_mask',#'slice_timings',
                                                   'discard', 'dof_mc', 'fwhm', 'cost_mc',
                                                   'bet_frac', 'robust', 'wm_thresh', 'dof_f', 'bbr_type', 
                                                   'interp', 'cost', 'bins', 'iso', 'bbr',
                                                   'T1w', 'mean_img', 'slice_corrected',
                                                   'mask', 'warplater', 'pipeline_names']
        
        
        #fields.append([key[0] for key in im_var if key != 'PIPELINE'])
        
        inputnode = Node(IdentityInterface(fields=fields, mandatory_inputs=False), name='inputnode_im')
        inputnode.iterables = im_var['ExtractROI'] + im_var['smoother'] + co_var['BET'] + co_var['registration']
        inputnode.synchronize = True
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['mc_mean_img', 'slice_time_corrected', 
                                                 'smoothed', 'outliers', 'plots', 'mc_par',
                                                 'warped_mean', 'warped', 'brain', 'out_mat', 'warp_file', 'intermediate_outputs', 'segmentations']), name='outnode', joinsource='inputnode_im', joinfield=['mc_mean_img', 'slice_time_corrected','smoothed', 'outliers', 'plots', 'mc_par', 'intermediate_outputs', 'segmentations'])
        else:
            outnode = Node(IdentityInterface(fields=['mc_mean_img', 'slice_time_corrected', 
                                                 'smoothed', 'outliers', 'plots', 'mc_par',
                                                 'warped_mean', 'warped', 'brain', 'out_mat', 'warp_file', 'intermediate_outputs', 'segmentations']), name='outnode')
        
        reg = self.coreg_flow()
        
        #skip dummy scans
        extract = Node(ExtractROI(t_size=-1, output_type='NIFTI_GZ'), name='extract')
        #motion correction (lots of other options that can be adjusted)
        #NOTE: altered file /opt/anaconda3/lib/python3.8/site-packages/nipype/interfaces/fsl/preprocess.py
        #      line 936 to add or LooseVersion(Info.version()) > LooseVersion("6.0.3") as it appears fsl 
        #      changed the file extension output for mean_vol option to match that of fsl 5
        
        #dof=dof_mc, cost=cost_mc, 
        
        mc = Node(MCFLIRT(mean_vol=True, save_plots=True, output_type='NIFTI_GZ'), name='mc')
        #slice timing correction
        slicetimer = Node(SliceTimer(index_dir=False, interleaved=True, output_type='NIFTI_GZ'), name='slicetimer')
        
        def smoother(base_dir, warped, susan, fwhm, frac_mask):
            #WORKFLOW ADAPTED FROM: https://nipype.readthedocs.io/en/latest/users/examples/fmri_fsl.html
            from nipype import Workflow, Node, IdentityInterface
            from nipype.interfaces.fsl import (BET, ImageMaths, ImageStats, SUSAN, IsotropicSmooth)
            import os, glob
            from functions import get_bright_thresh, getthreshop
            
            smooth = Workflow('smooth')
            smooth.base_dir = os.getcwd() #base_dir
            #smooth.config['execution']['remove_unnecessary_outputs'] = 'False'
            
            inputnode = Node(IdentityInterface(fields=['in_file']), name='inputnode')
            inputnode.inputs.in_file = warped
            
            outnode = Node(IdentityInterface(fields=['smoothed']), name='outnode')
            
            if susan:
                meanfunc = Node(ImageMaths(op_string='-Tmean', suffix='_mean'), name='meanfunc')
                meanfuncmask = Node(BET(mask=True, no_output=True, frac=frac_mask), name='meanfuncmask')
                maskfunc = Node(ImageMaths(suffix='_bet', op_string='-mas'), name='maskfunc')
                getthresh = Node(ImageStats(op_string='-p 2 -p 98'), name='getthreshold')
                threshold = Node(ImageMaths(out_data_type='char', suffix='_thresh'), name='threshold')
                
                
                medianval = Node(ImageStats(op_string='-k %s -p 50'), name='medianval')
                
                smooth_su = Node(SUSAN(fwhm=fwhm), name='smooth_su')
                
                smooth.connect([(inputnode, meanfunc, [('in_file', 'in_file')]),
                                   (meanfunc, meanfuncmask, [('out_file', 'in_file')]),
                                   (inputnode, maskfunc, [('in_file', 'in_file')]),
                                   (meanfuncmask, maskfunc, [('mask_file', 'in_file2')]),
                                   (maskfunc, getthresh, [('out_file', 'in_file')]),
                                   (maskfunc, threshold, [('out_file', 'in_file')]),
                                   (getthresh, threshold, [(('out_stat', getthreshop), 'op_string')]),
                                   (inputnode, medianval, [('in_file', 'in_file')]),
                                   (threshold, medianval, [('out_file', 'mask_file')]),
                                   (medianval, smooth_su, [(('out_stat', get_bright_thresh), 'brightness_threshold')]),
                                   (inputnode, smooth_su, [('in_file', 'in_file')]),
                                   (smooth_su, outnode, [('smoothed_file', 'smoothed')]),
                                ])
                
                smooth.run()
                smoothed = glob.glob(os.getcwd() + '/smooth/smooth_su/*.nii*')[0]
            else:
                smooth_iso = Node(IsotropicSmooth(fwhm=fwhm, output_type='NIFTI'), name='smooth_iso')
                smooth.connect(inputnode, 'in_file', smooth_iso, 'in_file')
                smooth.connect(smooth_iso, 'out_file', outnode, 'smoothed')
                smooth.run()
                
                smoothed = glob.glob(os.getcwd() + '/smooth/smooth_iso/*.nii*')[0]
                 
            return smoothed, glob.glob(os.getcwd() + '/smooth/**', recursive=True)
        
        smoothnode = Node(Function(input_names=['base_dir', 'warped', 'susan', 'fwhm', 'frac_mask'],
                                   output_names=['smooth', 'files'], function=smoother), name='smoothnode')
        smoothnode.inputs.base_dir = os.getcwd() #self.base_dir
        
        #NOT FSL NATIVE - may delete: Artifact Detection - determines outliers in functional images
        #PROBABLY REPLACE WITH MELODIC
        art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=3,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  plot_type='svg'),
                   name="art")
        
        imageproc.connect([(inputnode, reg, [('bet_frac', 'inputnode_reg.bet_frac'),
                                              ('robust', 'inputnode_reg.robust'),
                                              ('wm_thresh', 'inputnode_reg.wm_thresh'),
                                              ('dof_f', 'inputnode_reg.dof_f'),
                                              ('bbr_type', 'inputnode_reg.bbr_type'),
                                              ('interp', 'inputnode_reg.interp'),
                                              ('cost', 'inputnode_reg.cost'),
                                              ('bins', 'inputnode_reg.bins'),
                                              ('iso', 'inputnode_reg.iso'),
                                              ('bbr', 'inputnode_reg.bbr'),
                                              ('T1w', 'inputnode_reg.T1w'),
                                              ('warplater', 'inputnode_reg.warplater'),
                                              ('mask', 'inputnode_reg.mask')]),
                           (mc, reg, [('mean_img', 'inputnode_reg.mean_img')]),
                           (slicetimer, reg, [('slice_time_corrected_file', 'inputnode_reg.slice_corrected')]),
                           (reg, smoothnode, [('outnode.warped', 'warped')]),
                           (reg, outnode, [('outnode.warped_mean', 'warped_mean'),
                                            ('outnode.warped', 'warped'),
                                            ('outnode.brain', 'brain'),
                                            ('outnode.out_mat', 'out_mat'),
                                            ('outnode.warp_file', 'warp_file'),
                                            ('outnode.segmentations', 'segmentations')]),
                           ])
        
        
        imageproc.connect([(inputnode, slicetimer, [('TR', 'time_repetition')]),
                           (inputnode, extract, [('discard', 't_min')]),
                           (inputnode, extract, [('bold', 'in_file')]),
                           (reg, art, [('outnode.warped', 'realigned_files')]),
                           (extract, mc, [('roi_file', 'in_file')]),
                           (mc, slicetimer, [('out_file', 'in_file')]),
                           (mc, art, [('par_file', 'realignment_parameters')]),
                           
                           (inputnode, smoothnode, [('susan', 'susan'),
                                                    ('fwhm', 'fwhm'),
                                                    ('frac_mask', 'frac_mask')]),
                           
                           #(inputnode, slicetimer, [('slice_timings', 'custom_timings')]),
                           (slicetimer, outnode, [('slice_time_corrected_file', 'slice_time_corrected')]), #might pass in slice timing file
                           (mc, outnode, [('mean_img', 'mc_mean_img')]),
                           (mc, outnode, [('par_file', 'mc_par')]),
                           (art, outnode, [('outlier_files', 'outliers')]),
                           (art, outnode, [('plot_files', 'plots')]),
                           (smoothnode, outnode, [('smooth', 'smoothed'),
                                                  ('files', 'intermediate_outputs')]),
                           ])
        
        def out(base_dir, pipeline_st, task, st, sm, mc):
            from nipype.interfaces import DataSink
            from nipype import Node
            import re, os
            sink = Node(DataSink(base_directory=base_dir, parameterization=False), name='sink')
            folder_name = task
            session = re.search('/(_sessions_[A-Za-z0-9]+)/', sm)
            run = re.search('/(_runs_[0-9]+)/', sm)
            subject = re.search('/(_subject_[0-9A-Za-z]+/)', sm).group(1)
            pipeline = re.search('/mapflow/[A-Za-z_]+([0-9]+)', os.getcwd()).group(1)
            pipeline = 'pipeline_' + str(int(pipeline) + pipeline_st)
            
            sink.inputs.container = pipeline
            
            if session:
                folder_name += '/' + session.group(1) + '/'
            if run:
                folder_name += '/' + run.group(1) + '/'
            if folder_name == task:
                folder_name += '/'
                
            folder_name += 'preprocess/' + subject
            
            setattr(sink.inputs, folder_name + '.@st', st)
            setattr(sink.inputs, folder_name + '.@sm', sm)
            setattr(sink.inputs, folder_name + '.@mc', mc)
            
            sink.run()
        
        write = MapNode(Function(input_names=['base_dir', 'pipeline_st', 'task', 'st', 'sm', 'mc'],
                                 function=out), name='write', iterfield=['st', 'sm', 'mc'])
        write.inputs.base_dir = self.base_dir
        write.inputs.pipeline_st = self.pipeline
        write.inputs.task = self.task
        
        imageproc.connect([(outnode, write, [('slice_time_corrected', 'st'),
                                             ('mc_par', 'mc'),
                                             ('smoothed', 'sm')
                                             ]),
                           ])

        return imageproc
    
    #bet_frac=0.5, wm_thresh=0.5, dof_f=6, bbr=True, bbr_type='signed', interp='spline', iso=4, cost='mutualinfo', bins=640
    #robust=True
    def coreg_flow(self): #iso -> iso_resample
        #NOTE: this code is adapted from the nipype tutorial on preprocessing https://miykael.github.io/nipype_tutorial/notebooks/example_preprocessing.html 
        from functions import get_wm
        coregwf = Workflow(name='coregwf')
        coregwf.base_dir = os.getcwd() #self.base_dir
        
        fields=['bet_frac', 'robust', 'wm_thresh', 'dof_f', 'bbr_type', 
                                                   'interp', 'cost', 'bins', 'iso', 'bbr',
                                                   'T1w', 'mean_img', 'slice_corrected', 'mask', 'warplater', 'ind']
        
        
        #fields.append([key[0] for key in co_var if key != 'PIPELINE'])
        
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_reg')
        
        #for tup in inputs:
        #    setattr(inputnode.inputs, tup[0], tup[1])
            #in_name = 'inputnode.inputs.' + tup[0]
            #vars()[in_name] = tup[1]
        
        #inputnode.inputs = co_var['BET'] + co_var['registration']
        
        outnode = Node(IdentityInterface(fields=['warped_mean', 'warped', 'brain', 'out_mat', 'warp_file', 'segmentations', 'intermediate_outputs']), name='outnode')
         
        #Brain extraction (maybe frac is a parameter that can be adjusted, could alter -g vertical gradient intensity threshold)
        bet_anat = Node(BET(output_type='NIFTI_GZ'), name="bet_anat")#, iterfield='in_file')
        
        segment = Node(FAST(output_type='NIFTI_GZ'), name='segment')
        
        def registration(base_dir, T1w, mean_img, slice_corrected, bet, wm_file, bbr, wm_thresh, dof_f, bbr_type, interp, iso, cost, bins, warp_file, mask, warplater):
            from nipype.interfaces.fsl import (FAST, FLIRT, Threshold)
            from nipype import IdentityInterface
            from nipype import Workflow, Node, Function
            from functions import parse_xml
            from nipype.interfaces.fsl import ImageMeants, ConvertXFM, ApplyWarp, ExtractROI
            from os.path import join as opj
            import os, glob
            from functions import get_wm
            reg = Workflow(name='reg')
            reg.base_dir = os.getcwd()
            #reg.config['execution']['remove_unnecessary_outputs'] = 'False'
            #CHECK OPTIONS
            
            #MOVE SEGMENT OUTSIDE (SO INTERMEDIATES STORED), AND BBR SPECIFIC INTO BBR IF STATEMENT
            #segment = Node(FAST(in_files=bet, output_type='NIFTI_GZ'), name='segment')
            
            #CHECK OPTIONS
            threshold = Node(Threshold(thresh=wm_thresh, args='-bin', output_type='NIFTI_GZ'), name="threshold")
            threshold.inputs.in_file = wm_file
            #FLIRT HAS BINS AS INPUT #cost=cost,
            reg_pre = Node(FLIRT(dof=dof_f, in_file=mean_img, reference=bet, output_type='NIFTI_GZ'), name='reg_pre')
            #bbr_type=bbr_type,
            reg_bbr = Node(FLIRT(dof=dof_f, cost='bbr', reference=T1w, 
                                 in_file=mean_img, output_type='NIFTI_GZ',
                                 schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch')), name='reg_bbr')
            
            
            applywarp = Node(FLIRT(interp=interp, in_file=slice_corrected, reference=bet, apply_isoxfm=iso, output_type='NIFTI_GZ'), name='applywarp')
            applywarp_mean = Node(FLIRT(interp=interp, in_file=mean_img, reference=bet, apply_isoxfm=iso, output_type='NIFTI_GZ'), name='applywarp_mean')
            
            outnode =  Node(IdentityInterface(fields=['warped_mean', 'warped','out_mat']), name='outnode')
            
            if bbr:
                reg.connect([#(segment, threshold, [(('partial_volume_files', get_wm), 'in_file')]),
                             (threshold, reg_bbr, [('out_file', 'wm_seg')]),
                             (reg_pre, reg_bbr, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, applywarp, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, applywarp_mean, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, outnode, [('out_matrix_file', 'out_mat')]),
                             ])
                node_reg = 'reg_bbr'
            else:
                reg.connect([(reg_pre, applywarp, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_pre, applywarp_mean, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_pre, outnode, [('out_matrix_file', 'out_mat')]),
                             ])
                node_reg = 'reg_pre'
                
            if not warplater:
                to_std = Node(ApplyWarp(ref_file=mask, field_file=warp_file), name='to_std')
                reg.connect(applywarp, 'out_file', to_std, 'in_file')
                reg.connect(to_std, 'out_file', outnode, 'warped')
                node_warp = 'to_std'
            else:
                reg.connect(applywarp, 'out_file', outnode, 'warped')
                node_warp = 'applywarp'
                
            reg.connect(applywarp_mean, 'out_file', outnode, 'warped_mean')
            reg.run()
                
            out_mat = glob.glob(os.getcwd() + '/reg/' + node_reg + '/*.mat')[0]
            warped_mean = glob.glob(os.getcwd() + '/reg/applywarp_mean/*.nii*')[0]
            warped = glob.glob(os.getcwd() + '/reg/' + node_warp + '/*.nii*')[0]
            
            return out_mat, warped_mean, warped, glob.glob(os.getcwd() + '/reg/**', recursive=True)
        
        regnode = Node(Function(input_names=['base_dir', 'T1w', 'mean_img', 'slice_corrected', 'bet', 'wm_file', 'bbr', 
                                             'wm_thresh', 'dof_f', 'bbr_type', 'interp', 'iso', 'cost', 'bins',
                                             'warp_file', 'mask', 'warplater'],
                                output_names=['out_mat', 'warped_mean', 'warped', 'files'], 
                                function=registration), 
                                name='regnode')
        regnode.inputs.base_dir = os.getcwd() #self.base_dir
        
        warp_generator = self.get_warps()
        coregwf.connect(inputnode, 'mask', warp_generator, 'inputnode.ref_file')
        coregwf.connect(bet_anat, 'out_file', warp_generator, 'inputnode.brain')
        coregwf.connect(warp_generator, 'outnode.warp', regnode, 'warp_file')
        coregwf.connect([(warp_generator, outnode, [('outnode.warp', 'warp_file')]),
                         ])
        
        #coregwf.connect(inputnode, 'bbr', buff, 'bbr')
        #node connection
        coregwf.connect([(inputnode, regnode, [('bbr', 'bbr'),
                                               ('wm_thresh', 'wm_thresh'),
                                               ('dof_f', 'dof_f'),
                                               ('bbr_type', 'bbr_type'),
                                               ('interp', 'interp'),
                                               ('iso', 'iso'),
                                               ('cost', 'cost'),
                                               ('bins', 'bins'),
                                               ('mean_img', 'mean_img'),
                                               ('slice_corrected', 'slice_corrected'),
                                               #('warp_file', 'warp_file'),
                                               ('T1w', 'T1w'),
                                               ('warplater', 'warplater')]),
                         (inputnode, regnode, [('mask', 'mask')]),
                         #(inputnode, regnode, [('slice_corrected', 'slice_corrected')]),
                         #(inputnode, regnode, [('mean_img', 'applywarp_mean.in_file')]),
                         #(inputnode, regnode, [('mean_img', 'reg_pre.in_file')]),
                         (inputnode, bet_anat, [('bet_frac', 'frac')]),
                         (inputnode, bet_anat, [('robust', 'robust')]),
                         (inputnode, bet_anat, [('T1w', 'in_file')]),
                         #REGISTRATION
                         (bet_anat, regnode, [('out_file', 'bet')]),
                         #(bet_anat, regnode, [('out_file', 'applywarp.reference')]),
                         #(bet_anat, regnode, [('out_file', 'applywarp_mean.reference')]),
                         #(bet_anat, segment, [('out_file', 'in_files')]),
                         (bet_anat, segment, [('out_file', 'in_files')]),
                         (segment, regnode, [(('partial_volume_files', get_wm), 'wm_file')]),
                         #OUTPUT
                         (regnode, outnode, [('out_mat', 'out_mat')]),
                         (regnode, outnode, [('warped_mean', 'warped_mean')]),
                         (regnode, outnode, [('warped', 'warped')]),
                         (regnode, outnode, [('files', 'intermediate_outputs')]),
                         (bet_anat, outnode, [('out_file', 'brain')]),
                         (segment, outnode, [('partial_volume_files', 'segmentations')]),
                         ])
        
        return coregwf
    
class level1(spatial_normalization):
    def __init__(self, base_dir):
        super().__init__()
        self.base_dir = base_dir
        #self.pipeline = 0
        
    def construct(self, l1_var, norm):
        level1 = Workflow('level1')
        level1.base_dir = os.getcwd() #self.base_dir
        
        featnode = self.first_level_flow(l1_var)
        to_std = self.apply_warps(norm)
        #connode = self.generate_contrasts()
        
        inputnode = Node(IdentityInterface(fields=['TR', 'event_file', 'outliers', 'mc_par',
                                                   
                                                   'task', 'mask', 'warp', 'brain',
                                                   
                                                   'smoothed', 'segmentations']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['feat_dir', 'contrast_names', 'cope', 'varcope', 'bold']), name='outnode')
        
        level1.connect([(inputnode, to_std, [('warp', 'inputnode_app.warp_file'),
                                             ('mask', 'inputnode_app.ref_file')]),
                        (featnode, to_std, [('outnode.feat_dir', 'inputnode_app.feat_dir')]),
                        (to_std, outnode, [('outnode.cope', 'cope'),
                                           ('outnode.varcope', 'varcope'),
                                           ('outnode.bold', 'bold')]),
                        ])
        
        level1.connect([(inputnode, featnode, [('TR', 'inputnode_l1.TR'),
                                           #('discard', 'inputnode_l1.discard'),
                                           #('HP', 'inputnode_l1.HP'),
                                           #('thresh', 'inputnode_l1.thresh'),
                                           #('serial_cor', 'inputnode_l1.serial_cor'),
                                           #('base_switch', 'inputnode_l1.base_switch'),
                                           #('gamma', 'inputnode_l1.gamma'),
                                           #('dict_opt', 'inputnode_l1.dict_opt'),
                                           #('base_val', 'inputnode_l1.base_val'),
                                           ('event_file', 'inputnode_l1.event_file'),
                                           ('outliers', 'inputnode_l1.outliers'),
                                           ('mc_par', 'inputnode_l1.mc_par'),
                                           ('smoothed', 'inputnode_l1.smoothed'),
                                           ('task', 'inputnode_l1.task'),
                                           #('resting', 'inputnode_l1.resting'),
                                           ('mask', 'inputnode_l1.mask'),
                                           ('warp', 'inputnode_l1.warp'),
                                           ('segmentations', 'inputnode_l1.segmentations'),
                                           #('warp_post_feat', 'inputnode_l1.warp_post_feat'),
                                           ('brain', 'inputnode_l1.brain')]),
                        (featnode, outnode, [('outnode.feat_dir', 'feat_dir')]),
                        (featnode, outnode, [('outnode.contrast_names', 'contrast_names')]),
                        ])
        
        return level1
        
    def generate_contrasts(self):
        gencon = Workflow('gencon')
        gencon.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['session_info']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['contrast_names', 'contrasts']), name='outnode')
        
        def contrasts(session_info):
            contrasts = []
            identities = []
            for info in session_info:
                condition_names = []
                for task in info['cond']:
                    condition_names.append(task['name'])
                    
                num_tasks = len(condition_names)
                ind_con = []
                ind_ident = []
                group_con = []
                group_ident = []
                if num_tasks > 1:
                    for i, condition in enumerate(condition_names):
                        weights_specific = [0] * num_tasks
                        weights_specific[i] = 1
                        ind_con.append([condition, 'T', condition_names, weights_specific])
                        ind_ident.append(condition)
                        new_cond = condition + ' > others'
                        weights_specific = [-1/(num_tasks - 1)] * num_tasks
                        weights_specific[i] = 1
                        group_con.append([new_cond, 'T', condition_names, weights_specific])
                        group_ident.append(new_cond)
                   
                    contrasts.append(['average', 'T', condition_names, [1/num_tasks]*num_tasks])
                    identities.append('average')
                   
                    contrasts += ind_con
                    identities += ind_ident
                    contrasts += group_con
                    identities += group_ident
                   
                    #contrasts.append(['activation', 'F', ind_con])
                    #identities.append('activation')
                    #contrasts.append(['differences', 'F', group_con])
                    #identities.append('differences')
                else:
                    contrasts.append([condition_names[0], 'T', condition_names, [1]])
                    identities.append(condition_names[0])
                
            return identities, contrasts
        
        contrasts = Node(Function(input_names=['session_info'], output_names=['identities', 'contrasts'],
                                  function=contrasts), name='contrasts')
        
        gencon.connect(inputnode, 'session_info', contrasts, 'session_info')
        gencon.connect(contrasts, 'identities', outnode, 'contrast_names')
        gencon.connect(contrasts, 'contrasts', outnode, 'contrasts')
        
        return gencon
        
    #, TR, discard=4, HP=128, thresh=1000, serial_cor=True, base_switch=False, gamma=False, dict_opt='gammasigma', base_val=3
    def first_level_flow(self, l1_var):
        #Generates model
        #Outputs session_info
        #TO CONNECT: functional_runs, bids_event_file or event_files or subject_info, maybe outlier_files (art outlier files), maybe realignment parameters from MC
        
        l1_analysis = Workflow('l1_analysis')
        l1_analysis.base_dir = os.getcwd() #self.base_dir
        
        
        rng = len(l1_var['PIPELINE'][0][1])
        
        fields=['TR', 'discard', 'HP', 'thresh', 'serial_cor',
                                                   'base_switch', 'gamma', 'dict_opt', 'base_val',
                                                   'event_file', 'outliers', 'mc_par',
                                                   'smoothed', 'segmentations',
                                                   'task', 'resting', 'mask', 'warp', 'warp_post_feat', 'brain', 'ind']
        
        #fields.append([key[0] for key in l1_var['FLAMEO'] if key != 'PIPELINE'])
        
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_l1')
        inputnode.iterables = l1_var['ind']
        
        for tup in l1_var['session_info']:
            setattr(inputnode.inputs, tup[0], tup[1])
        for tup in l1_var['l1d']:
            setattr(inputnode.inputs, tup[0], tup[1])
        
        inputnode.synchronize = True
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['feat_dir', 'session_info', 'contrast_names', 'ev_files', 'intermediate_outputs']), name='outnode', joinsource='inputnode_l1', joinfield=['feat_dir', 'session_info', 'contrast_names', 'ev_files', 'intermediate_outputs'])
        else:    
            outnode = Node(IdentityInterface(fields=['feat_dir', 'session_info', 'contrast_names', 'ev_files', 'intermediate_outputs']), name='outnode')
        
        #MODIFY SO CAN CREATE SEED BASED OFF OF COORDINATES AND RADIUS
        #PASS IN CSF MASK AND GET SIGNAL FROM THERE TO ADD AS REGRESSOR
        #USE GLOBAL AVERAGE SIGNAL AS REGRESSOR AS WELL
        def session_info(event_file, outliers, mc_par, TR, HP, smoothed, task, resting, mask, warp, warp_post_feat, segmentations):
            from nipype import Node, Workflow, IdentityInterface
            from versatile import SpecifyModelVersatile
            from nipype.interfaces.fsl import ImageMeants, ExtractROI, ImageMaths, ApplyWarp, Threshold
            from nipype.interfaces.fsl.maths import MathsCommand
            from functions import parse_xml
            import os, glob, re
            import numpy as np
            
            model = Node(SpecifyModelVersatile(input_units='secs', parameter_source='FSL'), name='model')
            model.inputs.outlier_files = outliers
            model.inputs.realignment_parameters = mc_par
            model.inputs.time_repetition = TR
            model.inputs.high_pass_filter_cutoff = HP
            model.inputs.functional_runs = smoothed
                
            if event_file:
                model.inputs.bids_event_file = event_file
                session_info = model.run().outputs.session_info
            elif 'rest' in task:#nipype.algorithms.misc PickAtlas(although i've kind of implemented this, might be nice to add multiple ROIs -> get user input for options and index them in GA, pass list of list)
                if resting['atlas'] != 'ROI':
                    file, index, name = parse_xml(resting['atlas'], resting['seed'], mask)
                    
                    get_seed = Node(ExtractROI(in_file=file, t_min=int(index), t_size=1), name='get_seed')
                    roi = get_seed.run().outputs.roi_file
                    
                    thr_seed = Node(ImageMaths(op_string='-thr {seed_thr} -bin'.format(seed_thr=resting['seed_thr'])), name='thr_seed')
                    thr_seed.inputs.in_file = roi
                    thr_file = thr_seed.run().outputs.out_file
                    suffix = ''.join([word[0] for word in name.split()])
                else:
                    create_seed = Node(MathsCommand(in_file=mask), name='create_seed')
                    create_seed.inputs.args = '-mul 0 -add 1 -roi {x} 1 {y} 1 {z} 1 0 1'.format(x=resting['x'], y=resting['y'], z=resting['z'])
                    seed = create_seed.run().outputs.out_file
                    
                    make_sphere = Node(MathsCommand(in_file=seed), name='make_sphere')
                    make_sphere.inputs.args = '-kernel sphere {radius} -fmean'.format(radius=resting['radius'])
                    
                    thr_seed = Node(ImageMaths(op_string='-bin'), name='thr_seed')
                    thr_seed.inputs.in_file = make_sphere.run().outputs.out_file
                    thr_file = thr_seed.run().outputs.out_file
                    suffix = '_x{x}_y{y}_z{z}_r{r}'.format(x=resting['x'], y=resting['y'], z=resting['z'], r=resting['radius'])
                
                smooth_mni = smoothed
                
                if warp_post_feat:
                    mni = Node(ApplyWarp(ref_file=mask, field_file=warp), name='mni')
                    mni.inputs.in_file = smoothed
                    smooth_mni = mni.run().outputs.out_file
                    
                ev_name = re.search('task-([a-zA-Z]+)_', smooth_mni).group(1) + suffix
                    
                mean_ts = Node(ImageMeants(in_file=smooth_mni, out_file=ev_name), name='mean_ts')
                    
                mean_ts.inputs.mask = thr_file
                time_series = [mean_ts.run().outputs.out_file]
                
                model.inputs.event_files = time_series
                session_info = model.run().outputs.session_info
                #THIS IS MISSING A THRESHOLD FOR WHAT QUANTIFIES CSF
                if resting['CSF']:
                    ev_name = 'AvgCSFSig'
                    threshold = Node(Threshold(thresh=resting['CSF_thresh'], args='-bin'), name='threshold')
                    threshold.inputs.in_file = segmentations[0]
                    
                    mni_csf = Node(ApplyWarp(ref_file=mask, field_file=warp), name='mni')
                    mni_csf.inputs.in_file = threshold.run().outputs.out_file
                    csf_m = mni_csf.run().outputs.out_file
                    
                    mean_csf = Node(ImageMeants(in_file=smooth_mni, out_file=ev_name), name='mean_csf')
                    
                    mean_csf.inputs.mask = csf_m
                    time_series = mean_csf.run().outputs.out_file
                    
                    session_info[0]['regress'].append({'name': ev_name, 'val': np.loadtxt(time_series, dtype=float)})
                if resting['GLOBAL']:
                    ev_name = 'AvgGlobalSig'
                    mean_sig = Node(ImageMeants(in_file=smooth_mni, out_file=ev_name), name='mean_sig')
                    
                    mean_sig.inputs.mask = mask
                    time_series = mean_sig.run().outputs.out_file
                    session_info[0]['regress'].append({'name': ev_name, 'val': np.loadtxt(time_series, dtype=float)})
            else:
                print('ERROR')
                
            return session_info
                    
        
        modelspec = Node(Function(input_names=['event_file', 'outliers', 'mc_par', 'TR', 'HP', 'smoothed', 'task', 'resting', 'mask', 'warp', 'warp_post_feat', 'brain', 'segmentations'],
                                  output_names=['session_info'], function=session_info), name='modelspec')       
        
        #THIS MIGHT BREAK
        def correct_task_timing(session_info, TR, discard):
            for i, info in enumerate(session_info):
                for j, task in enumerate(info['cond']):
                    try:
                        for k, num in enumerate(task['onset']):
                            session_info[i]['cond'][j]['onset'][k] = num - (TR * discard)
                    except:
                        return session_info
                        
            return session_info
        
        
        correction = Node(Function(input_names=['session_info', 'TR', 'discard'], output_names=['session_info'], 
                                  function=correct_task_timing), name='correction')
        
        #Generate EV files
        #basis_func options: dgamma
        #base_opt options: derivs(bool, or: dictionary gamma(derivs, gammasigma, gammadelay)
        #TO CONNECT: modelspec session_info output
        #OPTIONAL: contrasts, orthogonalization
        #OUTPUTS: ev_files, fsf_files
        
        def specify_base(base_dir, gamma, dict_opt, base_val, base_switch, TR, serial_cor, session_info, contrasts):
            from nipype import Workflow, Node, IdentityInterface, DataSink
            from versatile import Level1DesignVersatile
            import os, glob
            l1design = Workflow('l1design')
            l1design.base_dir = os.getcwd() #base_dir
            #l1design.config['execution']['remove_unnecessary_outputs'] = 'False'
            #inputnode = Node(IdentityInterface(fields=['dict_opt', 'base_val', 'base_switch', 'TR', 
            #                                           'serial_cor', 'session_info', 'contrasts']), name='inputnode')
            
            outnode = Node(IdentityInterface(fields=['fsf_files', 'ev_files']), name='outnode')
            #THIS ONE CAN BE CHANGED SO INPUT BASE AS DICTIONARY -> WON'T NEED FUNCTION ASPECT
            
            #NO BASE FOR RESTING STATE
            if gamma:
                l1d = Node(Level1DesignVersatile(bases={'gamma': {dict_opt: base_val}}, 
                                        interscan_interval=TR, session_info=session_info, contrasts=contrasts,
                                        model_serial_correlations=serial_cor), name='l1d')
            else:
                l1d = Node(Level1DesignVersatile(bases={'dgamma':{'derivs':base_switch}}, 
                                        interscan_interval=TR, session_info=session_info, contrasts=contrasts,
                                        model_serial_correlations=serial_cor), name='l1d')
            
            l1design.connect([#(inputnode, l1d, [('TR', 'interscan_interval'),
                              #                  ('serial_cor', 'model_serial_correlations'),
                              #                  ('contrasts', 'contrasts'),
                              #                  ('session_info', 'session_info')]),
                              (l1d, outnode, [('fsf_files', 'fsf_files'),
                                              ('ev_files', 'ev_files')]),
                              ])
            
            l1design.run()
            
            ev_files = glob.glob(os.getcwd() + '/l1design/l1d/ev*.txt')
            fsf_files = glob.glob(os.getcwd() + '/l1design/l1d/*.fsf')[0]
            
            return ev_files, fsf_files, glob.glob(os.getcwd() + '/l1design/**', recursive=True) #l1design
        
        l1d = Node(Function(input_names=['base_dir', 'gamma', 'dict_opt', 'base_val', 'base_switch', 
                                         'TR', 'serial_cor', 'session_info', 'contrasts'],
                            output_names=['ev_files', 'fsf_files', 'files'], function=specify_base), name='l1d')
        l1d.inputs.base_dir = os.getcwd() #self.base_dir
        
        #Generate design files
        #INPUTS: ev_files, fsf_files
        #OUTPUTS: design_file, con_file
        #featmod = Node(FEATModel(output_type='NIFTI_GZ'), name='featmod')
        
        #Run FEAT
        feat = Node(FEAT(), name='feat')
        
        connode = self.generate_contrasts()
        
        l1_analysis.connect([(modelspec, connode, [('session_info', 'inputnode.session_info')]),
                             (connode, l1d, [('outnode.contrasts', 'inputnode.contrasts')]),
                             (connode, outnode, [('outnode.contrast_names', 'contrast_names')]),
                             ])
        
        def buffer(outliers, mc_par, HP, smoothed, resting, warp_post_feat, discard, gamma, dict_opt, 
                   base_val, base_switch, serial_cor, segmentations, i, rng):
            if rng > 1:
                return outliers[i], mc_par[i], HP[i], smoothed[i], resting[i], warp_post_feat[i], discard[i], gamma[i], dict_opt[i], base_val[i], base_switch[i], serial_cor[i], segmentations[i]
            else:
                return outliers, mc_par, HP, smoothed, resting, warp_post_feat, discard, gamma, dict_opt, base_val, base_switch, serial_cor, segmentations
        
        buff = Node(Function(input_names=['outliers', 'mc_par', 'HP', 'smoothed', 'resting', 'warp_post_feat', 'discard', 'gamma', 'dict_opt', 'base_val', 'base_switch', 'serial_cor', 'segmentations', 'i', 'rng'],
                             output_names=['outliers', 'mc_par', 'HP', 'smoothed', 'resting', 'warp_post_feat', 'discard', 'gamma', 'dict_opt', 'base_val', 'base_switch', 'serial_cor', 'segmentations'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        l1_analysis.connect([(inputnode, buff, [('outliers', 'outliers'),
                                                ('mc_par', 'mc_par'),
                                                ('HP', 'HP'),
                                                ('smoothed', 'smoothed'),
                                                ('resting', 'resting'),
                                                ('warp_post_feat', 'warp_post_feat'),
                                                #('brain', 'brain'),
                                                ('discard', 'discard'),
                                                ('gamma', 'gamma'),
                                                ('dict_opt', 'dict_opt'),
                                                ('base_val', 'base_val'),
                                                ('base_switch', 'base_switch'),
                                                ('serial_cor', 'serial_cor'),
                                                ('segmentations', 'segmentations'),
                                                ('ind', 'i')
                                                ]),
                             ])
        
        l1_analysis.connect([(inputnode, modelspec, [('event_file', 'event_file')]),
                    (buff, modelspec, [('outliers', 'outliers')]), #
                    (buff, modelspec, [('mc_par', 'mc_par')]), #
                    (inputnode, modelspec, [('TR', 'TR')]),
                    (buff, modelspec, [('HP', 'HP')]), #
                    (buff, modelspec, [('smoothed', 'smoothed'),
                                       ('resting', 'resting'),
                                       ('warp_post_feat', 'warp_post_feat'),
                                       ('segmentations', 'segmentations')]), #
                    (inputnode, modelspec, [('task', 'task'),
                                            #('resting', 'resting'), #
                                            ('mask', 'mask'),
                                            ('warp', 'warp')]),
                                            #('warp_post_feat', 'warp_post_feat')]), #
                                            #('brain', 'brain')]), #
                    (modelspec, correction, [('session_info', 'session_info')]),
                    (inputnode, correction, [('TR', 'TR')]),
                    (buff, correction, [('discard', 'discard')]), #
                    (buff, l1d, [('gamma', 'gamma')]), #
                    (buff, l1d, [('dict_opt', 'dict_opt')]), #
                    (buff, l1d, [('base_val', 'base_val')]), #
                    (buff, l1d, [('base_switch', 'base_switch')]), #
                    (inputnode, l1d, [('TR', 'TR')]),
                    (buff, l1d, [('serial_cor', 'serial_cor')]), #
                    (connode, l1d, [('outnode.contrasts', 'contrasts')]),
                    (correction, l1d, [('session_info', 'session_info')]),
                    #(modelspec, contrasts, [('session_info', 'session_info')]),
                    #(l1d, featmod, [('outnode.fsf_files', 'fsf_file')]),
                    #(l1d, featmod, [('outnode.ev_files', 'ev_files')]),
                    (l1d, feat, [('fsf_files', 'fsf_file')]),
                    (feat, outnode, [('feat_dir', 'feat_dir')]),
                    (modelspec, outnode, [('session_info', 'session_info')]),
                    (l1d, outnode, [('ev_files', 'ev_files'),
                                    ('files', 'intermediate_outputs')]),
                    #(featmod, feat, [('design_file', 'fsf_file')]),
                    ])
        
        
        #def film_command_string(brightness_threshold, mask_size, autocorr_estimate_only=True, autocorr_noestimate=False, 
        #                        fit_armodel=False, multitaper_product=10, smooth_autocorr=True, 
        #                        threshold=1000, tukey_window=10, use_pava=False):
        # POTENTIAL TO ADD, PEOPLE GENERALLY DON'T CHANGE THESE IT SEEMS
        # FILMGLS IS EQUIVALENT TO FEAT
        
        #ABSURD LIST OF OPTIONAL INPUTS
        #IN: input data file
        #OPT: design_file, brightness_thresh & mask_size (from SUSAN)
        #OUT: corrections, dof_file, param_estimates, sigmasquareds, thresholdac
        #film = Node(FILMGLS(threshold=thresh, autocorr_noestimate=not serial_cor, output_type='NIFTI_GZ'), name='film')
        
        return l1_analysis
    
# =============================================================================
# class spatial_normalization:
#     def __init__(self, base_dir):
#         self.base_dir = base_dir
#         self.dec = True
#         
#     def construct(self):
#         if self.dec:
#             self.dec = False
#             return self.get_warps()
#         else:
#             self.dec = True
#             return self.apply_warps()
#     
#     def get_warps(self):
#         from nipype.interfaces.fsl import ConvertXFM, InvWarp
#         genwarps = Workflow('genwarps')
#         genwarps.base_dir = os.getcwd() #self.base_dir
#         
#         inputnode = Node(IdentityInterface(fields=['brain', 'ref_file']), name='inputnode')
#         outnode = Node(IdentityInterface(fields=['warp', 'invwarp']), name='outnode')
#         
#         prelim = Node(FLIRT(dof=12, output_type='NIFTI_GZ'), name='prelim')
#         warp = Node(FNIRT(field_file=True), name='warp')
#         
#         warp_inv = Node(InvWarp(), name='warp_inv')
#         
#         genwarps.connect([(inputnode, prelim, [('brain', 'in_file')]),
#                           (inputnode, prelim, [('ref_file', 'reference')]),
#                           (inputnode, warp, [('brain', 'in_file')]),
#                           (inputnode, warp, [('ref_file', 'ref_file')]),
#                           (prelim, warp, [('out_matrix_file', 'affine_file')]),
#                           (warp, warp_inv, [('field_file', 'warp')]),
#                           (inputnode, warp_inv, [('brain', 'reference')]),
#                           (warp, outnode, [('field_file', 'warp')]),
#                           (warp_inv, outnode, [('inverse_warp', 'invwarp')]),
#                           ])
#         
#         return genwarps
#         
#     def apply_warps(self):
#         appwarps = Workflow('appwarps')
#         appwarps.base_dir = os.getcwd()
#         
#         inputnode = Node(IdentityInterface(fields=['feat_dir', 'warp_file', 'ref_file', 'needwarp']), name='inputnode')
#         outnode = Node(IdentityInterface(fields=['cope', 'varcope', 'bold']), name='outnode')
#         
#         templates = {'cope': 'stats/cope*.nii.gz',
#                      'varcope': 'stats/varcope*.nii.gz',
#                      'bold': 'filtered_func_data.nii.gz'}
#         
#         selectfiles = Node(SelectFiles(templates, sort_filelist=True), name='selectfiles')
#         
#         def identity(cope, varcope, bold, needwarp):
#             from nipype.interfaces.fsl import ImageMaths
#             from nipype import Node
#             
#             if not needwarp:
#                 clear = Node(ImageMaths(op_string='-mul 0 -bin'), name='clear')
#                 clear.inputs.in_file = cope
#                 clear.run()
#                 cope = clear.outputs.out_file
#                 varcope = cope
#                 bold = cope
#             
#             return cope, varcope, bold
#         
#         ident = Node(Function(input_names=['cope', 'varcope', 'bold', 'needwarp'],
#                               output_names=['cope', 'varcope', 'bold'], function=identity), name='ident')
#         
#         def ret_files(cope_orig, varcope_orig, bold_orig, cope_warp, varcope_warp, bold_warp, needwarp):
#             if not needwarp:
#                 return cope_orig, varcope_orig, bold_orig
#             else:
#                 return cope_warp, varcope_warp, bold_warp
#             
#         ret = Node(Function(input_names=['cope_orig', 'varcope_orig', 'bold_orig', 'cope_warp', 
#                                          'varcope_warp', 'bold_warp', 'needwarp'],
#                             output_names=['cope', 'varcope', 'bold'], function=ret_files), name='ret')
#         
#         
#         applywarp_c = MapNode(ApplyWarp(), name='applywarp_c', iterfield=['in_file'])
#         applywarp_v = MapNode(ApplyWarp(), name='applywarp_v', iterfield=['in_file'])
#         applywarp_bold = Node(ApplyWarp(), name='applywarp_bold')
#         
#         appwarps.connect([(inputnode, applywarp_c, [('ref_file', 'ref_file')]),
#                           (inputnode, applywarp_v, [('ref_file', 'ref_file')]),
#                           (inputnode, applywarp_bold, [('ref_file', 'ref_file')]),
#                           (inputnode, applywarp_c, [('warp_file', 'field_file')]),
#                           (inputnode, applywarp_v, [('warp_file', 'field_file')]),
#                           (inputnode, applywarp_bold, [('warp_file', 'field_file')]),
#                           (selectfiles, ident, [('cope', 'cope')]),
#                           (selectfiles, ident, [('varcope', 'varcope')]),
#                           (selectfiles, ident, [('bold', 'bold')]),
#                           (inputnode, ident, [('needwarp', 'needwarp')]),
#                           (ident, applywarp_c, [('cope', 'in_file')]),
#                           (ident, applywarp_v, [('varcope', 'in_file')]),
#                           (ident, applywarp_bold, [('bold', 'in_file')]),
#                           (applywarp_c, ret, [('out_file', 'cope_warp')]),
#                           (applywarp_v, ret, [('out_file', 'varcope_warp')]),
#                           (applywarp_bold, ret, [('out_file', 'bold_warp')]),
#                           (selectfiles, ret, [('cope', 'cope_orig')]),
#                           (selectfiles, ret, [('varcope', 'varcope_orig')]),
#                           (selectfiles, ret, [('bold', 'bold_orig')]),
#                           (ret, outnode, [('cope', 'cope')]),
#                           (ret, outnode, [('varcope', 'varcope')]),
#                           (ret, outnode, [('bold', 'bold')]),
#                           ])
#         
#         return appwarps
# =============================================================================
    
class level2:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        
    def construct(self, num_sub, l2_var, suffix=''):
        from math import comb
        l2analysis = Workflow('l2analysis' + suffix)
        l2analysis.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'mode', 'subjects', 'split_half', 'mask', 'max']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['groups', 'copes', 'varcopes', 'flameo_stats', 'zstats']), name='outnode')
        
        groups = self.group_construction()
        
# =============================================================================
#         com = comb(num_sub, int(num_sub / 2))
#         if not com % 2:
#             num_groups = com
#         else:
#             num_groups = int(com / 2)
# =============================================================================
        num_groups = num_sub #now the number of groups

        if suffix:
            num_groups = 1
        
        level2 = self.second_level_flow(num_groups, suffix, l2_var)
        
        l2analysis.connect([(inputnode, groups, [('subjects', 'inputnode.subjects'),
                                                 ('split_half', 'inputnode.split_half'),
                                                 ('max', 'inputnode.max')]),
                            (inputnode, level2, [('copes', 'inputnode_l2'+suffix+'.copes'),
                                                 ('varcopes', 'inputnode_l2'+suffix+'.varcopes'),
                                                 #('mode', 'inputnode_l2.mode'),
                                                 ('mask', 'inputnode_l2'+suffix+'.mask')]),
                            (groups, level2, [('outnode.groups', 'inputnode_l2'+suffix+'.groups')]),
                            (groups, outnode, [('outnode.groups', 'groups')]),
                            (level2, outnode, [('outnode.copes', 'copes'),
                                               ('outnode.var_copes', 'varcopes'),
                                               ('outnode.flameo_stats', 'flameo_stats'),
                                               ('outnode.zstats', 'zstats')]),
                            ])
        
        return l2analysis
        
    #subjects, split_half=False    
    def group_construction(self):
        import itertools
        from nipype import Node, JoinNode, MapNode, Workflow, SelectFiles, IdentityInterface, Function
        groups = Workflow('groups')
        groups.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['subjects', 'split_half', 'max']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['groups']), name='outnode')
        
        def make_groups(subjects, split_half, ref):
            from functions import random_combination
            from math import comb
            group_container = []
            if split_half:
                com = comb(len(subjects), int(len(subjects)/2))
                if len(subjects) % 2:
                    com /= 2
                if com > ref:
                    loops = ref
                else:
                    loops = com
                
                unique = set()
                
                while len(unique) < int(loops):
                    group = random_combination(subjects)
                    unique.add(group)
                
                group_container += [[list(tup[0]), list(tup[1])] for tup in list(unique)]
            else:
                if type(subjects) == list:
                    group_container.append([subjects])
                else:
                    group_container.append([[subjects]])
            #else:
                #PLACEHOLDER
            #    group_container.append([subjects, ['Single Group']])
            
            #NOTE: SELECTED RANDOM GROUPS WON'T BE CONSISTENT ACROSS GENERATIONS -> SET A SEED?
            #PROBABLY FINE HONESTLY
            return group_container
        
        group = Node(Function(input_names=['subjects', 'split_half', 'ref'],
                              output_names=['group_container'], function=make_groups), name='group')
        
        groups.connect([(inputnode, group, [('subjects', 'subjects'),
                                            ('split_half', 'split_half'),
                                            ('max', 'ref')]),
                        (group, outnode, [('group_container', 'groups')]),
                        ])
        
        return groups
        
    
    def second_level_flow(self, iternum, suffix, l2_var):
        from nipype import Rename, JoinNode
        l2 = Workflow(name='l2')
        l2.base_dir = os.getcwd()
        
        rng = len(l2_var['PIPELINE'][0][1])
        
        fields = ['groups', 'copes', 'varcopes', 'mask', 'ind']
        fields += [key[0] for key in l2_var['FLAMEO']]
        
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_l2'+suffix)
        inputnode.iterables = l2_var['ind']
        inputnode.synchronize = True
        
        for tup in l2_var['FLAMEO']:
            setattr(inputnode.inputs, tup[0], tup[1])
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['copes', 'var_copes', 'flameo_stats', 'zstats']), name='outnode', joinsource='inputnode_l2'+suffix, joinfield=['copes', 'var_copes', 'flameo_stats', 'zstats'])
        else:    
            outnode = Node(IdentityInterface(fields=['copes', 'var_copes', 'flameo_stats', 'zstats']), name='outnode')
        
        def construction(groups, copes, varcopes):
            import re
            contrasts = len(copes[0])
            merge_contain_c = []
            merge_contain_v = []
            num_copes = []
            
            for group in groups:
                for subset in group:
                    cont_c = []
                    cont_v = []
                    cont_g = []
                    for i in range(contrasts):
                        cont_c.append([cope[i] for cope in copes if re.search('_subject_([0-9S]+)', cope[i]).group(1) in subset])
                        cont_v.append([varcope[i] for varcope in varcopes if re.search('_subject_([0-9S]+)', varcope[i]).group(1) in subset])
                        cont_g.append(len(cont_v[i]))
                    
                    num_copes.append(cont_g)
                    merge_contain_c.append(cont_c)
                    merge_contain_v.append(cont_v)       
        
            return merge_contain_c, merge_contain_v, num_copes
        
        
        construct = Node(Function(input_names=['groups', 'copes', 'varcopes'], 
                                  output_names=['merge_contain_c', 'merge_contain_v', 'num_copes'],
                                  function=construction), name='construct')            
        
        l2model = MapNode(L2Model(), name='l2_model', iterfield=['num_copes'], nested=True)
        #run_mode, sigma_dofs, outlier_iter
        #should binarize this mask!!
        flameo = MapNode(FLAMEO(), name='flameo', 
                         iterfield=['cope_file', 'var_cope_file', 'design_file', 't_con_file', 
                                    'cov_split_file'], nested=True)
        
        def indexer(in_files, index):
            return in_files[index]        
        
        index_copes = Node(Function(input_names=['in_files', 'index'], 
                                    ouput_names=['out'], 
                                    function=indexer), name='index_copes')
        
        index_var = Node(Function(input_names=['in_files', 'index'], 
                                  ouput_names=['out'], 
                                  function=indexer), name='index_var')
        
        merge_copes =  MapNode(Merge(dimension='t'), name='merge_copes', 
                               iterfield=['in_files'])
        
        merge_var = MapNode(Merge(dimension='t'), name='merge_var', 
                            iterfield=['in_files'])
        
        if iternum == 1 or suffix:
            def simulate_merge(merged_file):
                return [merged_file]
            
            index_copes.inputs.index = 0
            index_var.inputs.index = 0
            
            cope_join = Node(Function(input_names=['merged_file'], output_names=['merged_file'], function=simulate_merge), name='cope_join') #Node(IdentityInterface(fields=['merged_file']), name='cope_join')
            var_join = Node(Function(input_names=['merged_file'], output_names=['merged_file'], function=simulate_merge), name='var_join')#Node(IdentityInterface(fields=['merged_file']), name='var_join')
            
            def find_con(files):
                import re
                contrasts = []
                for file in files[0]:
                    add = re.search('_flameo([0-9]+)', file).group(1)
                    contrasts.append(str(int(add) + 1))
                return [contrasts]
            
            re_cope = MapNode(Rename(format_string="cope%(contrast)s", keep_ext=True), name='re_cope', iterfield=['in_file', 'contrast'], nested=True)
            re_varcope = MapNode(Rename(format_string="varcope%(contrast)s", keep_ext=True), name='re_varcope', iterfield=['in_file', 'contrast'], nested=True)
            re_zstat = MapNode(Rename(format_string="zstat%(contrast)s", keep_ext=True), name='re_zstat', iterfield=['in_file', 'contrast'], nested=True)
            
            l2.connect([(flameo, re_cope, [('copes', 'in_file'),
                                           (('copes', find_con), 'contrast')]),
                        (flameo, re_varcope, [('var_copes', 'in_file'),
                                              (('var_copes', find_con), 'contrast')]),
                        (flameo, re_zstat, [('zstats', 'in_file'),
                                            (('zstats', find_con), 'contrast')]),
                        (re_cope, outnode, [('out_file', 'copes')]),
                        (re_varcope, outnode, [('out_file', 'var_copes')]),
                        (re_zstat, outnode, [('out_file', 'zstats')]),
                        ])
            
            def buffer(copes, varcopes, run_mode, i, rng):
                if rng > 1:
                    new_c = []
                    new_v = []
                    for j in range(len(copes)):
                        new_c.append(copes[j][i])
                        new_v.append(varcopes[j][i])
                    return new_c, new_v, run_mode[i]
                else:
                    return copes, varcopes, run_mode
        else:
            iternode_l2 = Node(IdentityInterface(fields=['subgroup']), name='iternode_l2')
            iternode_l2.iterables = [('subgroup', range(iternum*2))]
            cope_join = JoinNode(IdentityInterface(fields=['merged_file']), name='cope_join', joinsource='iternode_l2', joinfield='merged_file')
            var_join = JoinNode(IdentityInterface(fields=['merged_file']), name='var_join', joinsource='iternode_l2', joinfield='merged_file')
            
            l2.connect([(iternode_l2, index_copes, [('subgroup', 'index')]),
                        (iternode_l2, index_var, [('subgroup', 'index')]),
                        (flameo, outnode, [('copes', 'copes'),
                                           ('var_copes', 'var_copes'),
                                           ('zstats', 'zstats')]),
                        ])
            
            
            def buffer(copes, varcopes, run_mode, i, rng):
                if rng > 1:
                    new_c = []
                    new_v = []
                    for j in range(len(copes)):
                        new_c.append(copes[j][i][0])
                        new_v.append(varcopes[j][i][0])
                    return new_c, new_v, run_mode[i]
                else:
                    return copes, varcopes, run_mode
        
        buff = Node(Function(input_names=['copes', 'varcopes', 'run_mode', 'i', 'rng'],
                             output_names=['copes', 'varcopes', 'run_mode'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        l2.connect([(inputnode, construct, [('groups', 'groups')]),
                    (inputnode, buff, [('copes', 'copes')]),
                    (inputnode, buff, [('varcopes', 'varcopes')]),
                    (inputnode, buff, [('run_mode', 'run_mode')]),
                    (inputnode, buff, [('ind', 'i')]),
                    
                    (buff, construct, [('copes', 'copes')]),
                    (buff, construct, [('varcopes', 'varcopes')]),
                    
                    (construct, index_copes, [('merge_contain_c', 'in_files')]),
                    (construct, index_var, [('merge_contain_v', 'in_files')]),
                    #(iternode, index_copes, [('subgroup', 'index')]),
                    #(iternode, index_var, [('subgroup', 'index')]),
                    
                    (index_copes, merge_copes, [('out', 'in_files')]),
                    (index_var, merge_var, [('out', 'in_files')]),
                    
                    (merge_copes, cope_join, [('merged_file', 'merged_file')]),
                    (merge_var, var_join, [('merged_file', 'merged_file')]),
                    
                    (construct, l2model, [('num_copes', 'num_copes')]),
                    
                    (cope_join, flameo, [('merged_file', 'cope_file')]),
                    (var_join, flameo, [('merged_file', 'var_cope_file')]),
                    
                    (l2model, flameo, [('design_mat', 'design_file'),
                                       ('design_con', 't_con_file'),
                                       ('design_grp', 'cov_split_file')]),
                    (inputnode, flameo, [('mask', 'mask_file')]),
                    (buff, flameo, [('run_mode', 'run_mode')]),
                    (flameo, outnode, [#('copes', 'copes'),
                                       #('var_copes', 'var_copes'),
                                       #('zstats', 'zstats'),
                                       ('stats_dir', 'flameo_stats')]),
                    ])
        
        return l2
    
    
class split_half:
    def __init__(self, base_dir):
        self.base_dir = base_dir
    
    def construct(self, num_pipelines):
        split = Workflow('split')
        split.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['zstats', 'mask', 'groups', 
                                                   'preproc_bold', 'covariate_frame']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['score', 'R', 'R_lst', 'P', 'P_lst']), name='outnode')
        
        stat = self.stat_maps(num_pipelines)
        pred = self.prediction(num_pipelines)
        dist = self.distance(num_pipelines)
        
        split.connect([(inputnode, stat, [('zstats', 'inputnode_r.zstats'),
                                          ('mask', 'inputnode_r.mask')]),
                       (inputnode, pred, [('groups', 'inputnode_p.groups'),
                                          ('preproc_bold', 'inputnode_p.preproc_bold'),
                                          ('covariate_frame', 'inputnode_p.covariate_frame'),
                                          ('mask', 'inputnode_p.mask')]),
                       (pred, dist, [('outnode.pred_mean', 'inputnode_d.P')]),
                       (stat, dist, [('outnode.repro_mean', 'inputnode_d.R')]),
                       (pred, outnode, [('outnode.pred', 'P_lst'),
                                        ('outnode.pred_mean', 'P')]),
                       (stat, outnode, [('outnode.repro', 'R_lst'),
                                        ('outnode.repro_mean', 'R')]),
                       (dist, outnode, [('outnode.score', 'score')]),
                       ])
        
        return split
        
    def distance(self, num_pipelines):
        dist = Workflow('dist')
        dist.base_dir = os.getcwd() #self.base_dir
        
        rng = len(num_pipelines[0][1])
        
        fields = ['R', 'P', 'ind']
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_d')
        inputnode.iterables = num_pipelines
        inputnode.synchronize = True
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['score']), name='outnode', joinsource='inputnode_d', joinfield=['score'])
        else:
            outnode = Node(IdentityInterface(fields=['score']), name='outnode')
        
        def calc(self, R, P):
            import numpy as np
            perfect = np.array([1, 1])
            pipe = np.array([R, P])
            score = np.linalg.norm(perfect - pipe)
            return 1 - score
        
        calculation = Node(Function(input_names=['R', 'P'],
                                    output_names=['score'], function=calc), name='calculation')
        
        def buffer(R, P, i, rng):
            if rng > 1:
                return R[i], P[i]
            else:
                return R, P
        
        buff = Node(Function(input_names=['R', 'P', 'i', 'rng'],
                             output_names=['R', 'P'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        dist.connect([(inputnode, buff, [('R', 'R'),
                                         ('P', 'P'),
                                         ('ind', 'i')]),
                      (buff, calculation, [('R', 'R'),
                                                ('P', 'P')]),
                      (calculation, outnode, [('score', 'score')]),
                      ])
        
        return dist
        
    def stat_maps(self, num_pipelines):
        from nilearn.plotting import plot_img_comparison as plt
        from nilearn.input_data import NiftiMasker
        import numpy as np
        
        reproducibility = Workflow('reproducibility')
        reproducibility.base_dir = os.getcwd() #self.base_dir
        
        rng = len(num_pipelines[0][1])
        
        fields = ['zstats', 'mask', 'ind']
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_r')
        inputnode.iterables = num_pipelines
        inputnode.synchronize = True
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['repro', 'repro_mean']), name='outnode', joinsource='inputnode_r', joinfield=['repro', 'repro_mean'])
        else:
            outnode = Node(IdentityInterface(fields=['repro', 'repro_mean']), name='outnode')
        #FIGURE OUT HOW TO REPRESS ALL GRAPH OUTPUTS
        def compare(stats, mask):
            from nilearn.input_data import NiftiMasker
            from nilearn.image import math_img
            from nilearn.plotting import plot_img_comparison as plot
            import numpy as np
            import matplotlib.pyplot as plt
            plt.interactive(False)
            masker = NiftiMasker(math_img('img > 0', img=mask))
            masker.fit()
            num_groups = int(len(stats)/2)
            out_stats = []
            
            for i in range(num_groups):
                out = plot(stats[2*i], stats[2*i+1], masker, plot_hist=False)
                out_stats += out
            
            return out_stats, np.mean(out_stats)
        
        compare = Node(Function(input_names=['stats', 'mask'],
                                output_names=['repro', 'repro_mean'],
                                function=compare), name='compare')
        
        def buffer(zstats, i, rng):
            if rng > 1:
                return zstats[i]
            else:
                return zstats
        
        buff = Node(Function(input_names=['zstats', 'i', 'rng'],
                             output_names=['zstats'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        
        reproducibility.connect([(inputnode, compare, [('mask', 'mask')]),
                                 (inputnode, buff, [('zstats', 'zstats'),
                                                    ('ind', 'i')]),
                                 (buff, compare, [('zstats', 'stats')]),
                                 (compare, outnode, [('repro', 'repro'),
                                                     ('repro_mean', 'repro_mean')]),
                                 ])
        
        return reproducibility
    
    def prediction(self, num_pipelines):
        import numpy as np
        import nibabel as nib
        pred = Workflow('pred')
        pred.base_dir = os.getcwd() #self.base_dir
        
        rng = len(num_pipelines[0][1])
        
        fields = ['groups', 'preproc_bold', 'covariate_frame', 'mask', 'ind']
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_p')
        inputnode.iterables = num_pipelines
        inputnode.synchronize = True
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['pred', 'pred_mean']), name='outnode', joinsource='inputnode_p', joinfield=['pred', 'pred_mean'])
        else:
            outnode = Node(IdentityInterface(fields=['pred', 'pred_mean']), name='outnode')
        
        def predict(groups, covariate, bold, mask):
            import re
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.multioutput import MultiOutputClassifier
            import pandas as pd
            import numpy as np
            import nibabel as nib
            covariate = pd.read_table(covariate)
            covariates = covariate.set_index('participant_id')
            #covariates = covariate.columns
            num_cov = len(covariates.columns)
            #MAY NEED TO CHANGE SO THAT ADDITIONAL COLUMNS IN PARTICIPANTS ARE USED AS INPUTS AND ONLY DISEASE IS PREDICTED
            pred = []
            #TEST SPEED ON OLD FORMAT AS WELL
            mask_img = nib.load(mask)
            mask_ind = np.nonzero(mask_img.get_fdata())
            
            t_size = nib.load(bold[0][0]).shape[-1]
            
            size = (int(np.ceil(np.ravel(bold).shape[0])), np.shape(mask_ind)[1]*t_size)
            dat = np.zeros(size).astype(np.int16)
            subs = []
            ind = 0
            for file_cont in bold:
                for file in file_cont:
                    sub = re.search('_subject_([0-9]+)', file).group(1)
                    subs.append(sub)
                    dat[ind, :] = nib.load(file).get_fdata()[mask_ind].reshape(1, -1)
                    ind += 1
                    
            if num_cov:
                clf = RandomForestClassifier(random_state=2021, n_estimators=50)
                #if num_cov > 1:
                #    clf = MultiOutputClassifier(rfc)
                #else:
                #    clf = rfc
                for group in groups:
                    X1 = group[0]
                    X2 = group[1]
                    y1 = []
                    y2 = []
                    
                    ind1 = []
                    ind2 = []
                    
                    index = []
                    for subgroup in group:
                        for s in subgroup:
                            index = 'sub-' + s
                            group_ind = [pos for pos, val in enumerate(subs) if val == s]
                            if s in X1:
                                ind1 += group_ind
                                y1.append(covariates.loc[index].to_list()*len(group_ind))
                            elif s in X2:
                                ind2 += group_ind
                                y2.append(covariates.loc[index].to_list()*len(group_ind))
                            
                    Y1 = np.ravel(np.char.array(y1))
                    Y2 = np.ravel(np.char.array(y2))
                    
                    dat1 = dat[tuple(ind1),...]
                    dat2 = dat[tuple(ind2),...]
                    
                    clf = clf.fit(dat1, Y1)
                    pred.append(clf.score(dat2, Y2))
    
    
                    clf = clf.fit(dat2, Y2)
                    pred.append(clf.score(dat1, Y1))
            else:
                print("NO PARTICIPANTS FILE")
                
            return pred, np.mean(pred)
        
        def buffer(bold, i, rng):
            if rng > 1:
                new_b = []
                for j in range(len(bold)):
                    new_b.append(bold[j][i])
                return new_b
            else:
                return bold
        
        buff = Node(Function(input_names=['bold', 'i', 'rng'],
                             output_names=['bold'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        predict = Node(Function(input_names=['groups', 'covariate', 'bold', 'mask'],
                             output_names=['pred', 'pred_mean'], function=predict), name='predict')
        
        pred.connect([(inputnode, predict, [('groups', 'groups'),
                                            ('covariate_frame', 'covariate'),
                                            ('mask', 'mask')]),
                      (inputnode, buff, [('preproc_bold', 'bold'),
                                           ('ind', 'i')]),
                      (buff, predict, [('bold', 'bold')]),
                      (predict, outnode, [('pred', 'pred'),
                                          ('pred_mean', 'pred_mean')]),
                      ])
        
        return pred
                            
from nipype.interfaces.fsl import MultipleRegressDesign

class level3:#(level2):
    def __init__(self, base_dir):
        self.base_dir = base_dir
        #super().__init__(exp_dir, working_dir)
        
    def construct(self, l3_var):
        l3analysis = Workflow('l3analysis')
        l3analysis.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['covariates', 'copes', 'varcopes',
                                                   'mask', 'subjects']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['copes', 'var_copes', 'zstats', 
                                                 'flameo_stats']), name='outnode')
        
        test = self.setup()
        analysis = self.third_level_flow(l3_var)
        
        l3analysis.connect([(inputnode, test, [('covariates', 'inputnode.covariate'),
                                               ('subjects', 'inputnode.subjects')]),
                            (inputnode, analysis, [('copes', 'inputnode_l3.copes'),
                                                   ('varcopes', 'inputnode_l3.varcopes'),
                                                   ('mask', 'inputnode_l3.mask')]),
                            (test, analysis, [('outnode.regressors', 'inputnode_l3.regressors'),
                                              ('outnode.contrasts', 'inputnode_l3.contrasts'),
                                              ('outnode.group_ids', 'inputnode_l3.group_ids')]),
                            (analysis, outnode, [('outnode.copes', 'copes'),
                                                 ('outnode.var_copes', 'var_copes'),
                                                 ('outnode.zstats', 'zstats'),
                                                 ('outnode.flameo_stats', 'flameo_stats')]),
                            ])
        
        return l3analysis
        
    def setup(self):
        #CREATE PAIRED AND UNPAIRED T-TEST SETUP
        groups = Workflow('groups')
        groups.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['covariate', 'subjects']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['regressors', 'contrasts', 'group_ids']), name='outnode')
        
        #ONLY IMPLEMENTS UNPAIRED T-TEST AS OF NOW
        #ASK ERIN ABOUT DEMEANING, ORTHOGONALIZATION
        def t_test(covariate, subjects):
            import pandas as pd
            import numpy as np
            from functions import find_orth
            #PROBABLY PASS IN COVARIATE AS FILE NAME
            covariate = pd.read_table(covariate)
            covariates = covariate.set_index('participant_id')
            covariates = covariates.loc[['sub-' + sub for sub in subjects]]
            categories = covariates.columns
            groupcat = 'group'
            EVs = {}
            contrasts = []
            
            #IF NO GROUPS OR ANYTHING
            if len(categories) > 0:
                if groupcat not in categories:
                    groupcat = categories[0]
                    
                #ASSUMES unpaired t-test
                #THIS IS EXPECTING STRING CATEGORIES FOR GROUPS -> could probably eliminate need with inclusing of json file
                group = covariates.groupby(groupcat)
                num_groups = len(group.count())
                group_ids = (group.ngroup() + 1).to_list()
                encoded = group.ngroup().to_list()
                labels = covariates[groupcat].unique()
                contrast = []
                
                for i in range(num_groups):
                    ev = [1 if val == i else 0 for val in encoded]
                    EVs[labels[i]] = ev
                    
                    solo = (labels[i] + ' mean', 'T', [labels[i]], [1]) #-> I THINK THIS MIGHT BE AN F CONTRAST
                    contrast = [(labels[i] + '-' + lab, 'T', [labels[i], lab], [1,-1]) if lab != labels[i] else solo for lab in labels]
                    contrasts += contrast #contrasts.append(contrast)
                    
                #NOTE: FOR THE NON-GROUP COVARIATES THEY ARE ADDED AS IS RIGHT NOW -> NO DEMEANING/ORTHOGONALIZATION
                cov = covariates.drop(groupcat, axis=1)
                cat = categories.drop(groupcat)
                
                #orthog = np.array(list(EVs.values())).T
                EV_safe = EVs.copy()
                
                for c in cat:
                    labels = cov[c].unique()
                    if len(EV_safe) > 1:
                        for key in EV_safe:
                            if type(labels[0]) == str:
                                group = cov.groupby(c)
                                encoded = group.ngroup().to_list()
                                #orthog, vec = find_orth(orthog, np.array(encoded, ndmin=2).T)
                                #EVs[c] = vec.tolist()
                                EVs[c + '_' + key] = [encoded[i] if val else 0 for i, val in enumerate(EV_safe[key])]
                            else:
                                reg = cov[c].to_list()
                                EVs[c + '_' + key] = [reg[i] if val else 0 for i, val in enumerate(EV_safe[key])]
                                #orthog, vec = find_orth(orthog, np.array(cov[c].to_list(), ndmin=2).T)
                                #EVs[c] = vec.tolist()
            else:
                single_group = [1] * len(covariates.index)
                label = 'group_mean'
                EVs[label] = single_group
                group_ids = single_group
                contrasts.append([label, 'T', [label], [1]])
                
                    
            return EVs, contrasts, group_ids
        
        unpaired = Node(Function(input_names=['covariate', 'subjects'],
                                 output_names=['EVs', 'contrasts', 'group_ids'],
                                 function=t_test), name='unpaired')
        
        groups.connect([(inputnode, unpaired, [('covariate', 'covariate'),
                                               ('subjects', 'subjects')]),
                        (unpaired, outnode, [('EVs', 'regressors'),
                                             ('contrasts', 'contrasts'),
                                             ('group_ids', 'group_ids')]),
                        ])
        
        return groups

    def third_level_flow(self, l3_var):
        from nipype.interfaces.fsl import MultipleRegressDesign#, Randomise
        l3 = Workflow('l3')
        l3.base_dir = os.getcwd() #self.base_dir
        
        rng = len(l3_var['PIPELINE'][0][1])
        
        fields = ['regressors', 'contrasts', 'group_ids', 'copes', 'varcopes', 'mask', 'ind']#run_mode
        fields += [key[0] for key in l3_var['FLAMEO']]
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_l3')
        inputnode.iterables = l3_var['ind']
        inputnode.synchronize = True
        
        for tup in l3_var['FLAMEO']:
            setattr(inputnode.inputs, tup[0], tup[1])
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['copes', 'var_copes', 'zstats', 'flameo_stats']), name='outnode', joinsource='inputnode_l3', joinfield=['copes', 'var_copes', 'zstats', 'flameo_stats'])
        else:
            outnode = Node(IdentityInterface(fields=['copes', 'var_copes', 'zstats', 'flameo_stats']), name='outnode')
        
        model = Node(MultipleRegressDesign(), name='model')
        
        #NOTE: THIS MIGHT REQUIRE FURTHER SEPARATION INTO GROUPS INSTEAD OF JUST GROUPING BY CONTRASTS
        def group_contrast(lst):
            subs = len(lst)
            contrasts = len(lst[0])
            grouped = []
            for con in range(contrasts):
                grouped.append([lst[sub][con] for sub in range(subs)])
            
            return grouped
        
        merge_copes = MapNode(Merge(dimension='t'), name='merge_copes', iterfield=['in_files'])
        merge_var = MapNode(Merge(dimension='t'), name='merge_var', iterfield=['in_files'])
        
        flameo = MapNode(FLAMEO(), name='flameo', iterfield=['cope_file', 'var_cope_file'])
        
        #correct = MapNode(Randomise(), name='correct', iterfield=[])
        
        #NOTE: CAN MAKE MASK FILES USING LOWER LEVEL STATS MASK OUTPUT
        
        #TODO: FINISH CONNECTING TO FLAMEO, IMPLEMENT CLUSTER CORRECTION -> START BUILDING OUTER WORKFLOW
        
        def buffer(copes, varcopes, run_mode, i, rng):
            if rng > 1:
                new_c = []
                new_v = []
                for k in range(len(copes)):
                    new_c.append(copes[k][i][0])
                    new_v.append(varcopes[k][i][0])
                    
                return new_c, new_v, run_mode[i]
            else:
                return copes[0], varcopes[0], run_mode
        
        buff = Node(Function(input_names=['copes', 'varcopes', 'run_mode', 'i', 'rng'],
                             output_names=['copes', 'varcopes', 'run_mode'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        l3.connect([(inputnode, model, [('regressors', 'regressors'),
                                        ('contrasts', 'contrasts'),
                                        ('group_ids', 'groups')]),
                    (inputnode, buff, [('copes', 'copes'),
                                         ('varcopes', 'varcopes'),
                                         ('run_mode', 'run_mode'),
                                         ('ind', 'i')]),
                    (buff, merge_copes, [(('copes', group_contrast), 'in_files')]),
                    (buff, merge_var, [(('varcopes', group_contrast), 'in_files')]),
                    (inputnode, flameo, [('mask', 'mask_file')]),
                    (buff, flameo, [('run_mode', 'run_mode')]),
                    (merge_copes, flameo, [('merged_file', 'cope_file')]),
                    (merge_var, flameo, [('merged_file', 'var_cope_file')]),
                    (model, flameo, [('design_con', 't_con_file'),
                                     ('design_grp', 'cov_split_file'),
                                     ('design_fts', 'f_con_file'),
                                     ('design_mat', 'design_file')]),
                    (flameo, outnode, [('copes', 'copes'),
                                       ('var_copes', 'var_copes'),
                                       ('zstats', 'zstats'),
                                       ('stats_dir', 'flameo_stats')]),
                    ])
        
        return l3
    
from nipype.interfaces.base import CommandLine
from nipype.interfaces.fsl import Cluster, BinaryMaths, Threshold, SmoothEstimate
import re
#import subprocess import run
    
class correction:
    #ADAPTED FROM https://github.com/poldracklab/ds003-post-fMRIPrep-analysis/blob/master/workflows.py
    def __init__(self, base_dir):
        self.base_dir = base_dir
        
    def construct(self, corr):
        corr = corr['decision']
        corrected = Workflow('corrected')
        corrected.base_dir = os.getcwd() #self.base_dir
        
        rng = len(corr['PIPELINE'][0][1])
        
        fields = ['zstat', 'mask', 'copes', 'ind']
        fields += [key[0] for key in corr['DEC']]
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode')
        inputnode.iterables = corr['ind']
        inputnode.synchronize = True
        
        for tup in corr['DEC']:
            setattr(inputnode.inputs, tup[0], tup[1])
        
        outnode = Node(IdentityInterface(fields=['corrected']), name='outnode')
        
# =============================================================================
#         fwe = self.FWE()
#         clust = self.clusterFWE()
#         fdr = self.FDR()
# =============================================================================

        #WILL NEED TO DO SOME RENAMING OF THE FINAL OUTPUTS
        
        def decision(method, p, zstat, mask, connectivity, copes, z_thresh):
            import os, glob
            base_dir = os.getcwd()
            if method == 'fwe':
                from functions import FWE
                cor = FWE(base_dir, zstat, p, mask)
                cor.run()
                corrected = glob.glob(os.getcwd() + '/FWE/fwe_thresh/*.nii*')[0]
            elif method == 'cluster':
                from functions import clusterFWE
                cor = clusterFWE(base_dir, p, z_thresh, mask, connectivity, copes, zstat)
                cor.run()
                corrected = glob.glob(os.getcwd() + '/Cluster/cluster_all/*.nii*')[0]
            elif method == 'fdr':
                from functions import FDR
                cor = FDR(base_dir, zstat, p, mask)
                cor.run()
                corrected = glob.glob(os.getcwd() + '/FDR/corrected/*.nii*')[0]
            else:
                print('Cluster correction method not implemented')
                
            return corrected#, glob.glob(os.getcwd() + '/reg/**', recursive=True)
        
        def buffer(zstat, copes, method, connectivity, z_thresh, p, i, rng):
            if rng > 1:
                return zstat[i], copes[i], method[i], connectivity[i], z_thresh[i], p[i]
            else:
                return zstat, copes, method, connectivity, z_thresh, p
        
        buff = Node(Function(input_names=['zstat', 'copes', 'method', 'connectivity', 'z_thresh', 'p', 'i', 'rng'],
                             output_names=['zstat', 'copes', 'method', 'connectivity', 'z_thresh', 'p'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        
        dec = Node(Function(input_names=['method', 'p', 'zstat', 'mask', 'connectivity', 'copes', 'z_thresh'],
                            output_names='correction', function=decision), name='dec')
        #dec.iterables = corr
        

        corrected.connect([(inputnode, dec, [('mask', 'mask')]),
                           (inputnode, buff, [('zstat', 'zstat'),
                                             ('copes', 'copes'),
                                             ('method', 'method'),
                                             ('connectivity', 'connectivity'),
                                             ('z_thresh', 'z_thresh'),
                                             ('p', 'p'),
                                             ('ind', 'i')]),
                           (buff, dec, [('zstat', 'zstat'),
                                             ('copes', 'copes'),
                                             ('method', 'method'),
                                             ('connectivity', 'connectivity'),
                                             ('z_thresh', 'z_thresh'),
                                             ('p', 'p')]),
                           (dec, outnode, [('corrected', 'corrected')]),
                           ])
        
        return corrected
