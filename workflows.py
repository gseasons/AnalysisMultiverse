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

#mask_file = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain_mask.nii.gz')

#PREPROCESSING
class preprocess:
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        
    def construct(self):
        preprocess = Workflow('preprocess')
        preprocess.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['bold', 'T1w', 'susan', 'TR', 'frac_mask',#'slice_timings',
                                                   'discard', 'dof_mc', 'fwhm', 'cost_mc',
                                                   'bet_frac', 'robust', 'wm_thresh', 'dof_f', 'bbr_type', 
                                                   'interp', 'cost', 'bins', 'iso', 'bbr']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['smoothed', 'outliers', 'plots', 'mc_par',
                                                 'warped_mean', 'warped', 'brain', 'reg_out_mat']), name='outnode')
        
        im = self.improcess_flow()
        reg = self.coreg_flow()
        
        preprocess.connect([(inputnode, im, [('bold', 'inputnode.bold'),
                                             ('T1w', 'inputnode.T1w'),
                                             ('TR', 'inputnode.TR'),
                                             ('frac_mask', 'inputnode.frac_mask'),
                                             ('discard', 'inputnode.discard'),
                                             ('dof_mc', 'inputnode.dof_mc'),
                                             ('fwhm', 'inputnode.fwhm'),
                                             ('cost_mc', 'inputnode.cost_mc')]),
                            (inputnode, reg, [('bet_frac', 'inputnode.bet_frac'),
                                              ('robust', 'inputnode.robust'),
                                              ('wm_thresh', 'inputnode.wm_thresh'),
                                              ('dof_f', 'inputnode.dof_f'),
                                              ('bbr_type', 'inputnode.bbr_type'),
                                              ('interp', 'inputnode.interp'),
                                              ('cost', 'inputnode.cost'),
                                              ('bins', 'inputnode.bins'),
                                              ('iso', 'inputnode.iso'),
                                              ('bbr', 'inputnode.bbr'),
                                              ('T1w', 'inputnode.T1w')]),
                            
                            (im, reg, [('outnode.mc_mean_img', 'inputnode.mean_img'),
                                       ('outnode.slice_time_corrected', 'inputnode.slice_corrected')]),
                            (reg, im, [('outnode.warped', 'inputnode.warped')]),
                            
                            (im, outnode, [('outnode.smoothed', 'smoothed'),
                                           ('outnode.outliers', 'outliers'),
                                           ('outnode.plots', 'plots'),
                                           ('outnode.mc_par', 'mc_par')]),
                            (reg, outnode, [('outnode.warped_mean', 'warped_mean'),
                                            ('outnode.warped', 'warped'),
                                            ('outnode.brain', 'brain'),
                                            ('outnode.out_mat', 'reg_out_mat')]),
                            ])
        
        return preprocess
        
        
        #frac_mask=0.3, TR, options, susan=False, bbr=True, discard=4, dof_mc=6, fwhm=4, cost_mc='normcorr'
    def improcess_flow(self):
        #ASSUMES INTERLEAVED SLICE TIMING, USES NIPYPE ERROR DETECTION INSTEAD OF MELODIC
        imageproc = Workflow(name='imageproc')
        base_dir = opj(self.exp_dir, self.working_dir)
        imageproc.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['bold', 'T1w', 'warped', 'susan', 'TR', 'frac_mask',#'slice_timings',
                                                   'discard', 'dof_mc', 'fwhm', 'cost_mc']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['mc_mean_img', 'slice_time_corrected', 
                                                 'smoothed', 'outliers', 'plots', 'mc_par']), name='outnode')
        
        #skip dummy scans
        extract = Node(ExtractROI(t_size=-1, output_type='NIFTI'), name='extract')
        #motion correction (lots of other options that can be adjusted)
        #NOTE: altered file /opt/anaconda3/lib.python3.8/site-packages/nipype/interfaces/fsl/preprocess.py
        #      line 936 to add or LooseVersion(Info.version()) > LooseVersion("6.0.3") as it appears fsl 
        #      changed the file extension output for mean_vol option to match that of fsl 5
        
        #dof=dof_mc, cost=cost_mc, 
        
        mc = Node(MCFLIRT(mean_vol=True, save_plots=True, output_type='NIFTI'), name='mc')
        #slice timing correction
        slicetimer = Node(SliceTimer(index_dir=False, interleaved=True, output_type='NIFTI'), name='slicetimer')
        
        def smoother(base_dir, warped, susan, fwhm, frac_mask):
            #WORKFLOW ADAPTED FROM: https://nipype.readthedocs.io/en/latest/users/examples/fmri_fsl.html
            smooth = Workflow('smooth')
            smooth.base_dir = base_dir
            
            inputnode = Node(IdentityInterface(fields=['in_file']), name='inputnode')
            inputnode.inputs.in_file = warped
            
            outnode = Node(IdentityInterface(fields=['smoothed']), name='outnode')
            
            if susan:
                def get_bright_thresh(medianval):
                    return 0.75 * medianval
                
                meanfunc = Node(ImageMaths(op_string='-Tmean', suffix='_mean'), name='meanfunc')
                meanfuncmask = Node(BET(mask=True, no_output=True, frac=frac_mask), name='meanfuncmask')
                maskfunc = Node(ImageMaths(suffix='_bet', op_string='-mas'), name='maskfunc')
                getthresh = Node(ImageStats(op_string='-p 2 -p 98'), name='getthreshold')
                threshold = Node(ImageMaths(out_data_type='char', suffix='_thresh'), name='threshold')
                
                def getthreshop(thresh):
                    return '-thr %.10f -Tmin -bin' % (0.1 * thresh[1])
                
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
            else:
                 smooth_iso = Node(IsotropicSmooth(fwhm=fwhm, output_type='NIFTI'), name='smooth_iso')
                 smooth.connect(inputnode, 'in_file', smooth_iso, 'in_file')
                 smooth.connect(smooth_iso, 'out_file', outnode, 'smoothed')
                 
            return smooth
                 
        smoothnode = Node(Function(input_names=['base_dir', 'warped', 'susan', 'fwhm', 'frac_mask'],
                                   out_names=['smooth']), name='smoothnode')
        smoothnode.inputs.base_dir = base_dir
        
        #NOT FSL NATIVE - may delete: Artifact Detection - determines outliers in functional images
        #PROBABLY REPLACE WITH MELODIC
        art = Node(ArtifactDetect(norm_threshold=2,
                                  zintensity_threshold=3,
                                  mask_type='spm_global',
                                  parameter_source='FSL',
                                  use_differences=[True, False],
                                  plot_type='svg'),
                   name="art")
        
        imageproc.connect([(inputnode, slicetimer, [('TR', 'time_repetition')]),
                           (inputnode, extract, [('discard', 't_min')]),
                           (inputnode, extract, [('bold', 'in_file')]),
                           (inputnode, art, [('warped', 'realigned_files')]),
                           (extract, mc, [('roi_file', 'in_file')]),
                           (mc, slicetimer, [('out_file', 'in_file')]),
                           (mc, art, [('par_file', 'realignment_parameters')]),
                           
                           (inputnode, smoothnode, [('warped', 'warped'),
                                                    ('susan', 'susan'),
                                                    ('fwhm', 'fwhm'),
                                                    ('frac_mask', 'frac_mask')]),
                           
                           #(inputnode, slicetimer, [('slice_timings', 'custom_timings')]),
                           (slicetimer, outnode, [('slice_time_corrected_file', 'slice_time_corrected')]), #might pass in slice timing file
                           (mc, outnode, [('mean_img', 'mc_mean_img')]),
                           (mc, outnode, [('par_file', 'mc_par')]),
                           (art, outnode, [('outlier_files', 'outliers')]),
                           (art, outnode, [('plot_files', 'plots')]),
                           (smoothnode, outnode, [('outnode.smoothed', 'smoothed')]),
                           ])

        return imageproc
    
    
    
    
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
        
        outnode = Node(IdentityInterface(fields=['warped_mean', 'warped','brain', 'out_mat']), name='outnode')
        
        #Brain extraction (maybe frac is a parameter that can be adjusted, could alter -g vertical gradient intensity threshold)
        bet_anat = Node(BET(output_type='NIFTI_GZ'), name="bet_anat")#, iterfield='in_file')
        
        def registration(base_dir, T1w, mean_img, bbr, wm_thresh, dof_f, bbr_type, interp, iso, cost, bins):
            from nipype.interfaces.fsl import (FAST, FLIRT, Threshold)
            from nipype import Workflow, Node
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
    
#NOTE: SOME RECURSIVE NONSENSE IN CONNECTING FEATNODE AND CONNODE, MAY NOT WORK
class level1:
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        
    def construct(self):
        level1 = Workflow('level2')
        level1.base_dir = opj(self.exp_dir, self.working_dir)
        
        featnode = self.first_level_flow()
        connode = self.generate_contrasts()
        
        inputnode = Node(IdentityInterface(fields=['TR', 'discard', 'HP', 'thresh', 'serial_cor',
                                                   'base_switch', 'gamma', 'dict_opt', 'base_val',
                                                   'contrasts', 'event_file', 'outliers', 'mc_par',
                                                   'smoothed']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['feat_dir', 'contrast_names']), name='outnode')
        
        level1.connect([(inputnode, featnode, [('TR', 'inputnode.TR'),
                                           ('discard', 'inputnode.discard'),
                                           ('HP', 'inputnode.HP'),
                                           ('thresh', 'inputnode.thresh'),
                                           ('serial_cor', 'inputnode.serial_cor'),
                                           ('base_switch', 'inputnode.base_switch'),
                                           ('gamma', 'inputnode.gamma'),
                                           ('dict_opt', 'inputnode.dict_opt'),
                                           ('base_val', 'inputnode.base_val'),
                                           ('event_file', 'inputnode.event_file'),
                                           ('outliers', 'inputnode.outliers'),
                                           ('mc_par', 'inputnode.mc_par'),
                                           ('smoothed', 'inputnode.smoothed')]),
                        (featnode, connode, [('outnode.session_info', 'inputnode.session_info')]),
                        (connode, featnode, [('outnode.contrasts', 'inputnode.contrasts')]),
                        (featnode, outnode, [('outnode.feat_dir', 'feat_dir')]),
                        (connode, outnode, [('outnode.contrast_names', 'contrast_names')]),
                        ])
        
        return level1
        
    def generate_contrasts(self):
        gencon = Workflow('gencon')
        gencon.base_dir = opj(self.exp_dir, self.working_dir)
        
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
    def first_level_flow(self):
        #Generates model
        #Outputs session_info
        #TO CONNECT: functional_runs, bids_event_file or event_files or subject_info, maybe outlier_files (art outlier files), maybe realignment parameters from MC
        
        l1_analysis = Workflow('l1_analysis')
        base_dir = opj(self.exp_dir, self.working_dir)
        l1_analysis.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['TR', 'discard', 'HP', 'thresh', 'serial_cor',
                                                   'base_switch', 'gamma', 'dict_opt', 'base_val',
                                                   'contrasts', 'event_file', 'outliers', 'mc_par',
                                                   'smoothed']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['feat_dir', 'session_info']), name='outnode')
        
        modelspec = Node(SpecifyModel(input_units='secs', parameter_source='FSL'), name='modelspec')        
        
        def correct_task_timing(session_info, TR, discard):
            for i, info in enumerate(session_info):
                for j, task in enumerate(info['cond']):
                    for k, num in enumerate(task['onset']):
                        session_info[i]['cond'][j]['onset'][k] = num - (TR * discard)
                        
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
            l1design = Workflow('l1design')
            l1design.base_dir = base_dir
            
            inputnode = Node(IdentityInterface(fields=['dict_opt', 'base_val', 'base_switch', 'TR', 
                                                       'serial_cor', 'session_info', 'contrasts']), name='inputnode')
            
            outnode = Node(IdentityInterface(fields=['fsf_files', 'ev_files']), name='outnode')
            
            if gamma:
                l1d = Node(Level1Design(bases={'gamma': {dict_opt: base_val}}, 
                                        interscan_interval=TR, 
                                        model_serial_correlations=serial_cor), name='l1d')
            else:
                l1d = Node(Level1Design(bases={'dgamma':{'derivs':base_switch}}, 
                                        interscan_interval=TR, 
                                        model_serial_correlations=serial_cor), name='l1d')
            
            l1design.connect([(inputnode, l1d, [('TR', 'interscan_interval'),
                                                ('serial_cor', 'model_serial_correlations'),
                                                ('contrasts', 'contrasts'),
                                                ('session_info', 'session_info')]),
                              (l1d, outnode, [('fsf_files', 'fsf_files'),
                                              ('ev_files', 'ev_files')]),
                              ])
            
            return l1design
        
        l1d = Node(Function(input_names=['base_dir', 'gamma', 'dict_opt', 'base_val', 'base_switch', 
                                         'TR', 'serial_cor', 'session_info', 'contrasts'],
                            output_names=['l1design'], function=specify_base), name='l1d')
        l1d.inputs.base_dir = base_dir
        
        #Generate design files
        #INPUTS: ev_files, fsf_files
        #OUTPUTS: design_file, con_file
        #featmod = Node(FEATModel(output_type='NIFTI_GZ'), name='featmod')
        
        #Run FEAT
        feat = Node(FEAT(), name='feat')
        
        l1_analysis.connect([(inputnode, modelspec, [('event_file', 'bids_event_file')]),
                    (inputnode, modelspec, [('outliers', 'outlier_files')]),
                    (inputnode, modelspec, [('mc_par', 'realignment_parameters')]),
                    (inputnode, modelspec, [('TR', 'time_repetition')]),
                    (inputnode, modelspec, [('HP', 'high_pass_filter_cutoff')]),
                    (inputnode, modelspec, [('smoothed', 'functional_runs')]),
                    (modelspec, correction, [('session_info', 'session_info')]),
                    (inputnode, correction, [('TR', 'TR')]),
                    (inputnode, correction, [('discard', 'discard')]),
                    (inputnode, l1d, [('gamma', 'inputnode.gamma')]),
                    (inputnode, l1d, [('dict_opt', 'inputnode.dict_opt')]),
                    (inputnode, l1d, [('base_val', 'inputnode.base_val')]),
                    (inputnode, l1d, [('base_switch', 'inputnode.base_switch')]),
                    (inputnode, l1d, [('TR', 'inputnode.TR')]),
                    (inputnode, l1d, [('serial_cor', 'inputnode.serial_cor')]),
                    (inputnode, l1d, [('contrasts', 'inputnode.contrasts')]),
                    (correction, l1d, [('session_info', 'inputnode.session_info')]),
                    #(modelspec, contrasts, [('session_info', 'session_info')]),
                    #(l1d, featmod, [('outnode.fsf_files', 'fsf_file')]),
                    #(l1d, featmod, [('outnode.ev_files', 'ev_files')]),
                    (l1d, feat, [('outnode.fsf_files', 'fsf_file')]),
                    (feat, outnode, [('feat_dir', 'feat_dir')]),
                    (modelspec, outnode, [('session_info', 'session_info')]),
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
    
class spatial_normalization:
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        
    def construct(self):
        return self.spatial_normalization_flow()
        
    def spatial_normalization_flow(self):
        #from nipype.interfaces.fsl import FNIRT, ApplyWarp, FLIRT
        #from nipype import Node, MapNode, Workflow, SelectFiles, IdentityInterface
        #import os
        #from os.path import join as opj
        #COULD PROBABLY CHANGE THESE PARAMETERS - ASK ERIN, TONS FOR FNIRT
        
        warpflow = Workflow('warpflow')
        warpflow.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['feat_dir', 'brain', 'ref_file']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['cope', 'varcope', 'bold']), name='outnode')
        
        #ref_file = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
        #reference=ref_file
        prelim = Node(FLIRT(dof=12, output_type='NIFTI_GZ'), name='prelim')
        #WARP IS ACTING WEIRDLY, NOT OUTPUTTING WARP FILE
        warp = Node(FNIRT(field_file=True), name='warp')
        
        #applywarp_t = MapNode(ApplyWarp(ref_file=ref_file), name='applywarp_t', iterfield=['in_file'])
        #applywarp_z = MapNode(ApplyWarp(ref_file=ref_file), name='applywarp_z', iterfield=['in_file'])
        applywarp_c = MapNode(ApplyWarp(), name='applywarp_c', iterfield=['in_file'])
        applywarp_v = MapNode(ApplyWarp(), name='applywarp_v', iterfield=['in_file'])
        applywarp_bold = Node(ApplyWarp(), name='applywarp_bold')
        
        templates = {'cope': 'stats/cope*.nii.gz',
                     'varcope': 'stats/varcope*.nii.gz',
                     'bold': 'filtered_func_data.nii.gz'}
                     #'t_stat': 'stats/tstat*.nii.gz',
                     #'z_stat': 'stats/zstat*.nii.gz',
                     #}
        
        selectfiles = Node(SelectFiles(templates, sort_filelist=True), name='selectfiles')
        
        warpflow.connect([(inputnode, prelim, [('brain', 'in_file')]),
                          (inputnode, prelim, [('ref_file', 'reference')]),
                          (inputnode, warp, [('brain', 'in_file')]),
                          (inputnode, warp, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_c, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_v, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_bold, [('ref_file', 'ref_file')]),
                          (inputnode, selectfiles, [('feat_dir', 'base_directory')]),
                          (prelim, warp, [('out_matrix_file', 'affine_file')]),
                          #(warp, applywarp_t, [('field_file', 'field_file')]),
                          #(warp, applywarp_z, [('field_file', 'field_file')]),
                          (warp, applywarp_c, [('field_file', 'field_file')]),
                          (warp, applywarp_v, [('field_file', 'field_file')]),
                          #(selectfiles, applywarp_t, [('t_stat', 'in_file')]),
                          #(selectfiles, applywarp_z, [('z_stat', 'in_file')]),
                          (selectfiles, applywarp_c, [('cope', 'in_file')]),
                          (selectfiles, applywarp_v, [('varcope', 'in_file')]),
                          (selectfiles, applywarp_bold, [('bold', 'in_file')]),
                          (applywarp_c, outnode, [('out_file', 'cope')]),
                          (applywarp_v, outnode, [('out_file', 'varcope')]),
                          (applywarp_bold, outnode, [('out_file', 'bold')]),
                          ])
        
        return warpflow
    
class level2:
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        
    def construct(self, num_group):
        from math import comb
        l2analysis = Workflow('l2analysis')
        l2analysis.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['copes', 'varcopes', 'mode', 'subjects', 'split_half', 'mask']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['groups', 'copes', 'varcopes', 'flameo_stats', 'zstats']), name='outnode')
        
        groups = self.group_construction()
        
        com = comb(num_group, int(num_group / 2))
        if com % 2:
            num_groups = com
        else:
            num_groups = int(com / 2)
        
        level2 = self.second_level_flow(num_groups)
        
        l2analysis.connect([(inputnode, groups, [('subjects', 'inputnode.subjects'),
                                                 ('split_half', 'inputnode.split_half')]),
                            (inputnode, level2, [('copes', 'inputnode.copes'),
                                                 ('varcopes', 'inputnode.varcopes'),
                                                 ('mode', 'inputnode.mode'),
                                                 ('mask', 'inputnode.mask')]),
                            (groups, level2, [('outnode.groups', 'inputnode.groups')]),
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
        groups = Workflow('groups')
        base_dir = opj(self.exp_dir, self.working_dir)
        groups.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['subjects', 'split_half']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['groups']), name='outnode')
        
        def make_groups(base_dir, subjects, split_half):
            group_container = []
            if split_half:
                prelim = list(itertools.combinations(subjects, round(len(subjects)/2)))
                pre_len = len(prelim)
                for i in range(pre_len):
                    if not (pre_len % 2):
                        if i == (pre_len/2):
                            break
                        else:
                            group_container.append([list(prelim[i]), list(prelim[-(i+1)])])
                    else:
                        missed = [sub for sub in subjects if sub not in list(prelim[i])]
                        group_container.append([missed, list(prelim[i])])
            else:
                if type(subjects) == list:
                    group_container.append([subjects])
                else:
                    group_container.append([[subjects]])
            #else:
                #PLACEHOLDER
            #    group_container.append([subjects, ['Single Group']])
            
            return group_container
        
        group = Node(Function(input_names=['base_dir', 'subjects', 'split_half'],
                              out_names='group_container', function=make_groups), name='group')
        
        group.inputs.base_dir = base_dir
        
        groups.connect([(inputnode, group, [('subjects', 'subjects'),
                                                  ('split_half', 'split_half')]),
                        (group, outnode, [('group_container', 'groups')]),
                        ])
        
        return groups
        
    #layout.get(return_type='filename', extension='.tsv', suffix='participants')
    #import pandas as pd
    #B = pd.read_table(A[0])
    #B.columns
    
    #mode='fe'
    def second_level_flow(self, iternum):
        #from nipype import Node, JoinNode, MapNode, Workflow, SelectFiles, IdentityInterface, Function
        #from nipype.interfaces.fsl import L2Model, FLAMEO, Cluster, ImageMaths, Merge
        #import os
        #from os.path import join as opj
        
        #MAKE OUTNODE, PUT FLAMEO OUTSIDE OF ITERATION -> separate iteration from rest of inputnode
        #FLAMEO OUTPUTS -> cope, varcope files, zstat files (for input into split half -> create this next as well in different module)
        l2 = Workflow(name='l2')
        l2.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['groups', 'copes', 'varcopes', 'mode', 'mask']), name='inputnode')
        #FLAMEO STATS IS DIRECTORY WITH ALL OUTPUTS
        outnode = Node(IdentityInterface(fields=['copes', 'var_copes', 'flameo_stats', 'zstats']), name='outnode')
        
        iternode = Node(IdentityInterface(fields=['subgroup']), name='iternode')
        iternode.iterables = [('subgroup', range(iternum))]
        
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
                        cont_c.append([cope[i] for cope in copes if re.search('_subject_([0-9]+)', cope[i]).group(1) in subset])
                        cont_v.append([varcope[i] for varcope in varcopes if re.search('_subject_([0-9]+)', varcope[i]).group(1) in subset])
                        cont_g.append(len(subset))
                    
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
        
        cope_join = JoinNode(IdentityInterface(fields=['merged_file']), name='cope_join', joinsource='iternode', joinfield='merged_file')
        var_join = JoinNode(IdentityInterface(fields=['merged_file']), name='var_join', joinsource='iternode', joinfield='merged_file')
        
        l2.connect([(inputnode, construct, [('groups', 'groups')]),
                    (inputnode, construct, [('copes', 'copes')]),
                    (inputnode, construct, [('varcopes', 'varcopes')]),
                    
                    (construct, index_copes, [('merge_contain_c', 'in_files')]),
                    (construct, index_var, [('merge_contain_v', 'in_files')]),
                    (iternode, index_copes, [('subgroup', 'index')]),
                    (iternode, index_var, [('subgroup', 'index')]),
                    
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
                    (inputnode, flameo, [('mode', 'run_mode'),
                                         ('mask', 'mask_file')]),
                    (flameo, outnode, [('copes', 'copes'),
                                       ('var_copes', 'var_copes'),
                                       ('zstats', 'zstats'),
                                       ('stats_dir', 'flameo_stats')]),
                    ])
        
        return l2
    
    
class split_half:
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
    
    def construct(self):
        A = 3
        #'groups', 'preproc_bold',
        split = Workflow('split')
        split.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['zstats', 'mask', 'groups', 
                                                   'preproc_bold', 'covariate_frame']), name='inputnode')
        
        outnode = Node(IdentityInterface(fields=['score', 'R', 'R_lst', 'P', 'P_lst']), name='outnode')
        
        stat = self.stat_maps()
        pred = self.prediction()
        dist = self.distance()
        
        split.connect([(inputnode, stat, [('zstats', 'inputnode.zstats'),
                                          ('mask', 'inputnode.mask')]),
                       (inputnode, pred, [('groups', 'inputnode.groups'),
                                          ('preproc_bold', 'inputnode.preproc_bold'),
                                          ('covariate_frame', 'inputnode.covariate_frame')]),
                       (pred, dist, [('outnode.pred_mean', 'inputnode.P')]),
                       (stat, dist, [('outnode.repro_mean', 'inputnode.R')]),
                       (pred, outnode, [('outnode.pred', 'P_lst'),
                                        ('outnode.pred_mean', 'P')]),
                       (stat, outnode, [('outnode.repro', 'R_lst'),
                                        ('outnode.repro_mean', 'R')]),
                       (dist, outnode, [('outnode.score', 'score')]),
                       ])
        
        return split
        
    def distance(self):
        dist = Workflow('dist')
        dist.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['R', 'P']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['score']), name='outnode')
        
        def calc(self, R, P):
            import numpy as np
            perfect = np.array([1, 1])
            pipe = np.array([R, P])
            score = np.linalg.norm(perfect - pipe)
            return 1 - score
        
        calculation = Node(Function(input_names=['R', 'P'],
                                    output_names=['score'], function=calc), name='calculation')
        
        dist.connect([(inputnode, calculation, [('R', 'R'),
                                                ('P', 'P')]),
                      (calculation, outnode, [('score', 'score')]),
                      ])
        
        return dist
        
    def stat_maps(self):
        from nilearn.plotting import plot_img_comparison as plt
        from nilearn.input_data import NiftiMasker
        import numpy as np
        
        reproducibility = Workflow('reproducibility')
        reproducibility.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['zstats', 'mask']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['repro', 'repro_mean']), name='outnode')
        #FIGURE OUT HOW TO REPRESS ALL GRAPH OUTPUTS
        def compare(stats, mask):
            masker = NiftiMasker(mask)
            masker.fit()
            num_groups = int(len(stats)/2)
            out_stats = []
            for i in range(num_groups):
                out = plt(stats[2*i], stats[2*i+1], masker, plot_hist=False)
                out_stats += out
            
            return out_stats, np.mean(out_stats)
        
        compare = Node(Function(input_names=['stats', 'mask'],
                                output_names=['repro', 'repro_mean'],
                                function=compare), name='compare')
        
        reproducibility.connect([(inputnode, compare, [('zstats', 'stats'),
                                                       ('mask', 'mask')]),
                                 (compare, outnode, [('repro', 'repro'),
                                                     ('repro_mean', 'repro_mean')]),
                                 ])
        
        return reproducibility
    
    def prediction(self):
        import numpy as np
        import nibabel as nib
        pred = Workflow('pred')
        pred.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['groups', 'preproc_bold', 'covariate_frame']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['pred', 'pred_mean']), name='outnode')
        
        def predict(groups, covariate, bold):
            import re
            from sklearn.svm import LinearSVM
            from sklearn.multioutput import MultiOutputClassifier
            
            covariates = covariate.set_index('participant_id')
            #covariates = covariate.columns
            num_cov = len(covariates.columns)
            
            pred = []
            
            if num_cov:
                svm = LinearSVM(random_state=0)
                if num_cov > 1:
                    clf = MultiOutputClassifier(svm)
                else:
                    clf = svm
                    
                for group in groups:
                    X1 = group[0]
                    X2 = group[1]
                    X1_im = []
                    y1 = []
                    X2_im = []
                    y2 = []
                    for file in bold:
                        sub = re.search(file, '_subject_([0-9]+)').group(1)
                        index = 'sub-' + sub
                        if sub in X1:
                            img = np.array(nib.load(file)).reshape(1, -1)
                            X1_im.append(img)
                            y1.append(covariates.loc[index].to_list())
                        elif sub in X2:
                            img = np.array(nib.load(file)).reshape(1, -1)
                            X2_im.append(img)
                            y2.append(covariates.loc[index].to_list())
                        else:
                            print("ERROR")
                            
                    Y1 = np.char.array(y1)
                    Y2 = np.char.array(y2)
                    
                    dat1 = np.array(X1_im)
                    dat2 = np.array(X2_im)
                    
                    clf = clf.fit(dat1, Y1)
                    pred.append(clf.score(dat2, Y2))
                    
                    clf = clf.fit(dat2, Y2)
                    pred.append(clf.score(dat1, Y1))
            else:
                print("NO PARTICIPANTS FILE")
                
            return pred, np.mean(pred)
        
        predict = Node(Function(input_names=['groups', 'covariate', 'bold'],
                             output_names=['pred', 'pred_mean'], function=predict), name='predict')
        
        pred.connect([(inputnode, predict, [('groups', 'groups'),
                                            ('covariate_frame', 'covariate'),
                                            ('preproc_bold', 'bold')]),
                      (predict, outnode, [('pred', 'pred'),
                                          ('pred_mean', 'pred_mean')]),
                      ])
        
        return pred
                            
from nipype.interfaces.fsl import MultipleRegressDesign

class level3:#(level2):
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        #super().__init__(exp_dir, working_dir)
        
    def construct(self):
        l3analysis = Workflow('l3analysis')
        l3analysis.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['covariates', 'copes', 'varcopes',
                                                   'mask', 'mode']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['copes', 'var_copes', 'zstats', 
                                                 'flameo_stats']), name='outnode')
        
        test = self.setup()
        analysis = self.third_level_flow()
        
        l3analysis.connect([(inputnode, test, [('covariates', 'inputnode.covariate')]),
                            (inputnode, analysis, [('copes', 'inputnode.copes'),
                                                   ('varcopes', 'inputnode.varcopes'),
                                                   ('mask', 'inputnode.mask'),
                                                   ('mode', 'inputnode.mode')]),
                            (test, analysis, [('outnode.regressors', 'inputnode.regressors'),
                                              ('outnode.contrasts', 'inputnode.contrasts'),
                                              ('outnode.group_ids', 'inputnode.group_ids')]),
                            (analysis, outnode, [('outnode.copes', 'copes'),
                                                 ('outnode.var_copes', 'var_copes'),
                                                 ('outnode.zstats', 'zstats'),
                                                 ('outnode.flameo_stats', 'flameo_stats')]),
                            ])
        
        return l3analysis
        
    def setup(self):
        #CREATE PAIRED AND UNPAIRED T-TEST SETUP
        groups = Workflow('groups')
        groups.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['covariate']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['regressors', 'contrasts', 'group_ids']), name='outnode')
        
        #ONLY IMPLEMENTS UNPAIRED T-TEST AS OF NOW
        #ASK ERIN ABOUT DEMEANING, ORTHOGONALIZATION
        def t_test(self, covariate):
            covariates = covariate.set_index('participant_id')
            categories = covariate.columns
            groupcat = 'groups'
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
                    
                    solo = [labels[i] + ' mean', 'T', [labels[i]], [1]] #-> I THINK THIS MIGHT BE AN F CONTRAST
                    contrast = [[labels[i] + '-' + lab, 'T', [labels[i], lab], [1,-1]] if lab != labels[i] else solo for lab in labels]
                    contrasts.append(contrast)
                    
                #NOTE: FOR THE NON-GROUP COVARIATES THEY ARE ADDED AS IS RIGHT NOW -> NO DEMEANING/ORTHOGONALIZATION
                cov = covariates.drop(groupcat, axis=1)
                cat = categories.remove(groupcat)
                
                for c in cat:
                    labels = cov[c].unique()
                    
                    if type(labels[0]) == str:
                        encode = labels.ngroup().to_list()
                        EVs[c] = encode
                    else:
                        EVs[c] = cov[c].to_list()
            else:
                single_group = [1] * len(covariates.index)
                label = 'group_mean'
                EVs[label] = single_group
                group_ids = single_group
                contrasts.append([label, 'T', [label], [1]])
                
                    
            return EVs, contrasts, group_ids
        
        unpaired = Node(Function(input_names=['covariate'],
                                 output_names=['EVs', 'contrasts'],
                                 function=t_test), name='unpaired')
        
        groups.connect([(inputnode, unpaired, [('covariate', 'covariate')]),
                        (unpaired, outnode, [('EVs', 'regressors'),
                                             ('contrasts', 'contrasts'),
                                             ('group_ids', 'group_ids')]),
                        ])
        
        return groups

    def third_level_flow(self):
        from nipype.interfaces.fsl import MultipleRegressDesign#, Randomise
        l3 = Workflow('l3')
        l3.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['regressors', 'contrasts', 'group_ids',
                                                   'copes', 'varcopes', 'mask', 'mode']), name='inputnode')
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
        
        l3.connect([(inputnode, model, [('regressors', 'regressors'),
                                        ('contrasts', 'contrasts'),
                                        ('group_ids', 'groups')]),
                    (inputnode, merge_copes, [(('copes', group_contrast), 'in_files')]),
                    (inputnode, merge_var, [(('varcopes', group_contrast), 'in_files')]),
                    (inputnode, flameo, [('mask', 'mask_file'),
                                         ('mode', 'run_mode')]),
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
    def __init__(self, exp_dir, working_dir):
        self.exp_dir = exp_dir
        self.working_dir = working_dir
        
    def construct(self):
        A=3
    
    def FWE(self):
        A=3
        #TODO: NON CLUSTER FWE, GENERATE ANALYSIS WORKFLOW
        
    def clusterFWE(self):
        Cluster = Workflow('Cluster')
        Cluster.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['zstat', 'mask', 'connectivity', 'copes' 
                                                   'threshold', 'pthreshold']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['corrected']), name='outnode')
        
        smoothness = MapNode(SmoothEstimate(), name='smoothness')
        
        cluster_pos = MapNode(Cluster(out_index_file=True, out_local_max_txt_file=True), name='cluster_pos')
        cluster_neg = MapNode(Cluster(out_index_file=True, out_local_max_txt_file=True), name='cluster_neg')
        
        zstat_inv = MapNode(BinaryMaths(operation='mul', operand_value=-1), name='zstat_inv')
        cluster_inv = MapNode(BinaryMaths(operation='mul', operand_value=-1), name='cluster_inv')
        cluster_all = MapNode(BinaryMaths(operation='add'), name='cluster_all')
        
        Cluster.connect([(inputnode, smoothness, [('zstat', 'zstat_file'),
                                                  ('mask', 'mask')]),
                         (inputnode, cluster_neg, [('connectivity', 'connectivity'),
                                                   ('threshold', 'threshold'),
                                                   ('pthreshold', 'pthreshold'),
                                                   ('copes', 'cope_file')]),
                         (inputnode, cluster_pos, [('connectivity', 'connectivity'),
                                                   ('threshold', 'threshold'),
                                                   ('pthreshold', 'pthreshold'),
                                                   ('zstat', 'in_file'),
                                                   ('copes', 'cope_file')]),
                         (inputnode, zstat_inv, [('zstat', 'in_file')]),
                         (smoothness, cluster_neg, [('volume', 'volume'),
                                                    ('dlh', 'dlh')]),
                         (smoothness, cluster_pos, [('volume', 'volume'),
                                                    ('dlh', 'dlh')]),
                         (zstat_inv, cluster_neg, [('out_file', 'in_file')]),
                         (cluster_neg, cluster_inv, [('threshold_file', 'in_file')]),
                         (cluster_pos, cluster_all, [('threshold_file', 'in_file')]),
                         (cluster_inv, cluster_all, [('out_file', 'operand_file')]),
                         (cluster_all, outnode, [('out_file', 'corrected')]),
                         ])
        
        return Cluster
        
        
    def FDR(self):
        FDR = Workflow('FDR')
        FDR.base_dir = opj(self.exp_dir, self.working_dir)
        
        inputnode = Node(IdentityInterface(fields=['zstat', 'mask', 'p']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['corrected']), name='outnode')
        
        p_file = MapNode(ImageMaths(op_string='-ztop', suffix='_pval'), name='p_file', iterfield=['in_file'])
        
        def fdr(p_im, mask, q): #q = 0.05
            cmd = ('fdr -i {p_im} -m {mask} -q {q}')
            cl = CommandLine(cmd.format(p_im=p_im, mask=mask, q=q))
            results = cl.run().runtime.stdout
            thresh = 1 - float(re.search('([0-9\.]+)', results).group(1))
            
            form = '-mul -1 -add 1 -thr {thresh} -mas {mask}'.format(thresh=thresh, mask=mask)
            
            return form
        
        form = MapNode(Function(input_names=['p_im', 'mask', 'q'],
                               output_names=['form_str'], function=fdr), name='form', iterfield=['p_im'])
        
        
        corrected = MapNode(ImageMaths(suffix='_fdr'), name='corrected', iterfield=['in_file'])
        
        FDR.connect([(inputnode, p_file, [('zstat', 'in_file')]),
                     (inputnode, form, [('mask', 'mask'),
                                        ('p', 'q')]),
                     (p_file, form, [('out_file', 'p_im')]),
                     (form, corrected, [('form_str', 'op_string')]),
                     (p_file, corrected, [('out_file', 'in_file')]),
                     (corrected, outnode, [('out_file', 'corrected')])
                     ])
        
        return FDR
            
    
    #Randomize?
                
                
                
        
        

#A = preprocess('/Users/grahamseasons/fMRI/output_comp', 'working_dir').construct()
#B = level1('/Users/grahamseasons/fMRI/output_comp', 'working_dir').construct()
#C = spatial_normalization('/Users/grahamseasons/fMRI/output_comp', 'working_dir').construct()
#D = level2('/Users/grahamseasons/fMRI/output_comp', 'working_dir', 4).construct()
#E = split_half('/Users/grahamseasons/fMRI/output_comp', 'working_dir').construct()
F = level3('/Users/grahamseasons/fMRI/output_comp', 'working_dir').construct()
B=3