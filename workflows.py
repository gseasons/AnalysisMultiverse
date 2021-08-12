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


#PREPROCESSING
class preprocess:
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        
    def construct(self):
        A = 3
        
    def preprocess_flow(self, TR, options, susan=False, bbr=True, discard=4, dof_mc=6, fwhm=4, cost_mc='normcorr'):
        from os.path import join as opj
        import os
        import json
        from nipype.interfaces.fsl import (BET, ExtractROI, FAST, FLIRT, ImageMaths, ImageStats,
                                           MCFLIRT, SliceTimer, Threshold, SUSAN, IsotropicSmooth)
        #from nipype.interfaces.spm import Smooth
        from nipype.interfaces.utility import IdentityInterface
        from nipype.interfaces.io import SelectFiles, DataSink
        from nipype.algorithms.rapidart import ArtifactDetect
        from nipype import Workflow, Node
        
        from nipype.interfaces.spm import Smooth
        import nipype.interfaces.matlab as matlab
        from nipype.interfaces import spm
        spm.SPMCommand.set_mlab_paths(matlab_cmd='/Applications/MATLAB_R2017b.app/bin/matlab')
        
        matlab.MatlabCommand.set_default_matlab_cmd('/Applications/MATLAB_R2017b.app/bin/matlab')
        
        #MAIN PREPROCESSING WORKFLOW
        wf = self.coreg_flow(options)
        
        #skip dummy scans
        extract = Node(ExtractROI(t_min=discard, t_size=-1, output_type='NIFTI'), name='extract')
        #motion correction (lots of other options that can be adjusted)
        #NOTE: altered file /opt/anaconda3/lib.python3.8/site-packages/nipype/interfaces/fsl/preprocess.py
        #      line 936 to add or LooseVersion(Info.version()) > LooseVersion("6.0.3") as it appears fsl 
        #      changed the file extension output for mean_vol option to match that of fsl 5
        
        #dof=dof_mc, cost=cost_mc, 
        
        mc = Node(MCFLIRT(mean_vol=True, save_plots=True, output_type='NIFTI'), name='mc')
        #slice timing correction
        slicetimer = Node(SliceTimer(index_dir=False, interleaved=True, output_type='NIFTI', time_repetition=TR), name='slicetimer')
        #image smoothing - will need to figure out bright thresh later, might be 0.75 * median of run
        if susan:
            #WORKFLOW ADAPTED FROM: https://nipype.readthedocs.io/en/latest/users/examples/fmri_fsl.html
            def get_bright_thresh(medianval):
                return 0.75 * medianval
            
            io = Node(IdentityInterface(fields=['in_file']), name='io')
            io_o = Node(IdentityInterface(fields=['out_file']), name='io_o')
            meanfunc = Node(ImageMaths(op_string='-Tmean', suffix='_mean'), name='meanfunc')
            meanfuncmask = Node(BET(mask=True, no_output=True, frac=0.3), name='meanfuncmask')
            maskfunc = Node(ImageMaths(suffix='_bet', op_string='-mas'), name='maskfunc')
            getthresh = Node(ImageStats(op_string='-p 2 -p 98'), name='getthreshold')
            threshold = Node(ImageMaths(out_data_type='char', suffix='_thresh'), name='threshold')
            
            def getthreshop(thresh):
                print(thresh)
                return '-thr %.10f -Tmin -bin' % (0.1 * thresh[1])
            
            medianval = Node(ImageStats(op_string='-k %s -p 50'), name='medianval')
            
            def getbrightthresh(medianvals):
                return 0.75 * medianvals
            
            smooth_su = Node(SUSAN(fwhm=fwhm), name='smooth_su')
            
            smooth = Workflow(name='smooth')
            smooth.base_dir = opj(self.exp_dir, self.working_dir)
            smooth.connect([(io, meanfunc, [('in_file', 'in_file')]),
                               (meanfunc, meanfuncmask, [('out_file', 'in_file')]),
                               (io, maskfunc, [('in_file', 'in_file')]),
                               (meanfuncmask, maskfunc, [('mask_file', 'in_file2')]),
                               (maskfunc, getthresh, [('out_file', 'in_file')]),
                               (maskfunc, threshold, [('out_file', 'in_file')]),
                               (getthresh, threshold, [(('out_stat', getthreshop), 'op_string')]),
                               (io, medianval, [('in_file', 'in_file')]),
                               (threshold, medianval, [('out_file', 'mask_file')]),
                               (medianval, smooth_su, [(('out_stat', get_bright_thresh), 'brightness_threshold')]),
                               (io, smooth_su, [('in_file', 'in_file')]),
                               (smooth_su, io_o, [('smoothed_file', 'out_file')]),
                            ])
        
        else:
            smooth = Node(IsotropicSmooth(fwhm=fwhm, output_type='NIFTI'), name='smooth')#SUSAN(brightness_threshold=bright_thresh, fwhm=fwhm), name="smooth")
            #smooth = Node(Smooth(fwhm=fwhm), name='smooth')
        
        #NOT FSL NATIVE - may delete: Artifact Detection - determines outliers in functional images
        #PROBABLY REPLACE WITH MELODIC
        art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=3,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  plot_type='svg'),
                   name="art")
        
        preproc = Workflow(name='preproc')
        preproc.base_dir = opj(self.exp_dir, self.working_dir)
        
        preproc.connect([(extract, mc, [('roi_file', 'in_file')]),
                 (mc, slicetimer, [('out_file', 'in_file')]),
                 
                 (mc, wf, [('mean_img', 'reg_pre.in_file')]), #mean_img
                 (mc, wf, [('mean_img', 'applywarp_mean.in_file')]), #mean_img
                 
                 (slicetimer, wf, [('slice_time_corrected_file', 'applywarp.in_file')]), #might pass in slice timing file
                 
                 (wf, art, [('applywarp.out_file', 'realigned_files')]),
                 (mc, art, [('par_file', 'realignment_parameters')]),
                 ])
        
        if bbr:
            preproc.connect([(mc, wf, [('mean_img', 'reg_bbr.in_file')])]) #mean_img
            
        if susan:
            preproc.connect(wf, 'applywarp.out_file', smooth, 'io.in_file')
        else:
            preproc.connect(wf, 'applywarp.out_file', smooth, 'in_file')
            
        return preproc
    
    
    
    
    #bet_frac=0.5, wm_thresh=0.5, dof_f=6, bbr=True, bbr_type='signed', interp='spline', iso=4, cost='mutualinfo', bins=640
    #robust=True
    def coreg_flow(self): #iso -> iso_resample
        #NOTE: this code is adapted from the nipype tutorial on preprocessing https://miykael.github.io/nipype_tutorial/notebooks/example_preprocessing.html 
        coregwf = Workflow(name='coregwf')
        base_dir = opj(self.exp_dir, self.working_dir)
        coregwf.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['bet_frac', 'robust', 'wm_thresh', 'dof_f', 'bbr_type', 
                                                   'interp', 'cost', 'bins', 'iso', 'bbr',
                                                   'T1w', 'mean_img', 'slice_corrected']), name='inputnode')
        
        #Brain extraction (maybe frac is a parameter that can be adjusted, could alter -g vertical gradient intensity threshold)
        bet_anat = Node(BET(output_type='NIFTI_GZ'), name="bet_anat")#, iterfield='in_file')
        
        def registration(base_dir, T1w, mean_img, bbr, wm_thresh, dof_f, bbr_type, interp, iso, cost, bins):
            from nipype.interfaces.fsl import (FAST, FLIRT, Threshold)
            from nipype import Workflow, Node, MapNode
            reg = Workflow('reg')
            reg = base_dir
            
            #CHECK OPTIONS
            segment = Node(FAST(output_type='NIFTI_GZ'), name='segment')
            
            def get_wm(files):
                return files[-1]
            
            #CHECK OPTIONS
            threshold = Node(Threshold(thresh=wm_thresh, args='-bin', output_type='NIFTI_GZ'), name="threshold")
            
            #FLIRT HAS BINS AS INPUT
            reg_pre = Node(FLIRT(dof=dof_f, cost=cost, output_type='NIFTI_GZ'), name='reg_pre')
            
            reg_bbr = Node(FLIRT(dof=dof_f, cost='bbr', bbr_type=bbr_type, reference=T1w, 
                                 in_file=mean_img, output_type='NIFTI_GZ',
                                 schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch')), name='reg_bbr')
            
            
            applywarp = Node(FLIRT(interp=interp, apply_isoxfm=iso, output_type='NIFTI_GZ'), name='applywarp')
            applywarp_mean = Node(FLIRT(interp=interp, apply_isoxfm=iso, output_type='NIFTI_GZ'), name='applywarp_mean')
            
            outnode =  Node(IdentityInterface(fields=['warped_mean', 'warped','out_mat']), name='outnode')
            
            if bbr:
                reg.connect([(bet_anat, segment, [('out_file', 'in_files')]),
                             (segment, threshold, [(('partial_volume_files', get_wm), 'in_file')]),
                             (threshold, reg_bbr, [('out_file', 'wm_seg')]),
                             (reg_pre, reg_bbr, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, applywarp, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, applywarp_mean, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, outnode, [('out_matrix_file', 'out_mat')]),
                             ])
            else:
                reg.connect([(reg_pre, applywarp, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_pre, applywarp_mean, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_pre, outnode, [('out_matrix_file', 'out_mat')]),
                             ])
            
            reg.connect(applywarp_mean, 'out_file', outnode, 'warped_mean')
            reg.connect(applywarp, 'out_file', outnode, 'warped')
            
            return reg
        
        regnode = Node(Function(input_names=['base_dir', 'bbr', 'wm_thresh', 'dof_f', 
                                             'bbr_type', 'interp', 'iso', 'cost', 'bins',
                                             'T1w', 'mean_img'],
                                output_names=['reg'], 
                                function=registration), 
                                name='regnode')
        regnode.inputs.base_dir = base_dir
        
        
        outnode = Node(IdentityInterface(fields=['warped_mean', 'warped','brain', 'out_mat']), name='outnode')
                
    
        #node connection
        coregwf.connect([(inputnode, regnode, [('bbr', 'bbr'),
                                               ('wm_thresh', 'wm_thresh'),
                                               ('dof_f', 'dof_f'),
                                               ('bbr_type', 'bbr_type'),
                                               ('interp', 'interp'),
                                               ('iso', 'iso'),
                                               ('cost', 'cost'),
                                               ('bins', 'bins'),
                                               ('T1w', 'T1w'),
                                               ('mean_img', 'mean_img')]),
                         (inputnode, regnode, [('slice_corrected', 'applywarp.in_file')]),
                         (inputnode, regnode, [('mean_img', 'applywarp_mean.in_file')]),
                         (inputnode, regnode, [('mean_img', 'reg_pre.in_file')]),
                         (inputnode, bet_anat, [('bet_frac', 'frac')]),
                         (inputnode, bet_anat, [('robust', 'robust')]),
                         #REGISTRATION
                         (bet_anat, regnode, [('out_file', 'reg_pre.reference')]),
                         (bet_anat, regnode, [('out_file', 'applywarp.reference')]),
                         (bet_anat, regnode, [('out_file', 'applywarp_mean.reference')]),
                         #OUTPUT
                         (regnode, outnode, [('outnode.out_mat', 'out_mat')]),
                         (regnode, outnode, [('outnode.warped_mean', 'warped_mean')]),
                         (regnode, outnode, [('outnode.warped', 'warped')]),
                         (bet_anat, outnode, [('out_file', 'brain')]),
                         ])
        
        return coregwf