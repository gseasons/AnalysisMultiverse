#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 12:27:29 2021

@author: grahamseasons
"""

#TODO: OPTION TO JUST USE DEFAULT FMRI PREP PREPROCESSING
#DOES NOT HAVE SUPPORT FOR MULTIPLE RUNS
#WORKFLOW CONSTRUCTION VERY SLOW -> WORK TO IMPROVE
class fmri_multiverse:
    def __init__(self, data_path="/Users/grahamseasons/fMRI/test_data/ds000114/", saved_analysis=False):
        self.exp_dir = '/Users/grahamseasons/fMRI/output_comp'
        self.out_dir = 'datasink'
        self.working_dir = 'working_dir'
        self.data_path = data_path
        self.saved_analysis = saved_analysis
        self.get_data(data_path, saved_analysis)
        
    def get_data(self, path, analysis):
        from os.path import join as opj
        from bids.layout import BIDSLayout
        from nipype.interfaces.io import BIDSDataGrabber
        from nipype import Workflow, Node, JoinNode, Function, DataSink
        #GET AND RETURN IMPORTANT DATA, maybe write into containers pyradigm for organization
        #JUST DO ONE TASK FOR NOW, and only one session
        
        layout = BIDSLayout(self.data_path)
        subj = layout.get_subjects()
        types = layout.get_datatypes()
        session = layout.get_sessions()
        #im_type = layout.get_suffix(datatype='func') #change modality depending on info for anat, dwi, func
        tasks = layout.get_tasks()
        
        #CHANGE SO THAT RETURNS SLICE TIMINGS AS WELL
        #COULD TURN THIS INTO A NODE AS WELL -> A GOOD IDEA! TODO
        def metadata(task, ses, layout):
            #from bids.layout import BIDSLayout
            #layout = BIDSLayout("/Users/grahamseasons/fMRI/test_data/ds000114/")
            meta_path = layout.get(task=task, session=ses, datatype='func', suffix='bold', return_type='filename')[0]
            meta_data = layout.get_metadata(meta_path)
            TR = meta_data['RepetitionTime']
            return TR
        
        def remove_container(cont):
            if 'T1w' in cont[0] and 'sub-02' in cont[0]:
                return cont[0].replace('.gz', '')
            return cont[0]
        
        def event_grabber(file, path):
            import re
            from bids.layout import BIDSLayout
            layout = BIDSLayout(path)
            
            print('event grabber\n\n\n')
            print(file)
            
            task = re.search('task-([0-9A-Za-z]+)_bold', file).group(1)
            event_file = layout.get(task=task, extension='.tsv')
            if len(event_file) > 1:
                sub = re.search('/sub-([0-9]+)/', file).group(1)
                ses = re.search('_ses-([A-Za-z]+)_task', file).group(1)
                event_file = layout.get(task=task, session=ses, subject=sub, extension='.tsv')
            
            print(event_file[0])
                
            return event_file[0]
            
        event_node = Node(Function(input_names=['file', 'path'],
                                   output_names=['event_file'],
                                   function=event_grabber), name='event_node')
        
        event_node.inputs.path = self.data_path
        
            
        bbr = True
        susan=False
        pipelineID = "A"
            
        #greg = Node(Function(function=insight, input_names=["files"], output_names=[]), name='greg')
        #greg2 = Node(Function(function=met, input_names=["files"], output_names=[]), name='greg2')
        
        datasink = Node(DataSink(base_directory=self.exp_dir, container=self.out_dir), name="datasink")
        
        #substitutions = [('.nii_mean_reg','_mean_reg')]
        #datasink.inputs.substitutions = substitutions
        tasks = ['fingerfootlips']
        session = ['test']
            
        for task in tasks:
            for ses in session:
                TR = metadata(task, ses, layout)
                pre_wf = self.preprocess_flow(TR, options=3)
                
                if 'anat' and 'func' in types:
                    bids_dg = Node(BIDSDataGrabber(), name='bids_dg')
                    bids_dg.inputs.base_dir = self.data_path
                    bids_dg.inputs.output_query = {'bold_files': dict(datatype='func', suffix='bold', session=ses, task=task),
                                              #'bold_events': dict(datatype='func', suffix='events', session=ses),
                                              'T1w_files': dict(datatype='anat', suffix='T1w', session=ses)}
                    
                    bids_dg.iterables = ('subject', subj[:2])
                    
                    #WILL WANT TO PUT THESE IN WORKFLOW ASSEMBLY
                    data_proc = Workflow(name="data_proc")
                    data_proc.base_dir = opj(self.exp_dir, self.working_dir)
                    data_proc.connect([(bids_dg, pre_wf,[(('bold_files', remove_container), 'extract.in_file')]),
                                      (bids_dg, pre_wf, [(('T1w_files', remove_container), 'coregwf.bet_anat.in_file')]),
                                      
                                      (pre_wf, datasink, [('mc.par_file', 'preproc.' + pipelineID + '.@par')]),
                                      #smooth.smoothed_file
                                      #(pre_wf, datasink, [('smooth.out_file', 'preproc.' + pipelineID + '.@smooth')]),
                                      (pre_wf, datasink, [('coregwf.applywarp_mean.out_file', 'preproc.' + pipelineID + '.@mean')]),
                                      
                                      (pre_wf, datasink, [('coregwf.bet_anat.out_file', 'preproc.' + pipelineID + '.@brain')]),
                                      #OUTLIER FILES HOLD VOLUMES THAT HAVE ARTIFACTS, CAN READ THEM IN
                                      (pre_wf, datasink, [('art.outlier_files', 'preproc.' + pipelineID + '@outlier_files'),
                                                       ('art.plot_files', 'preproc.' + pipelineID + '@plot_files')]),
                                      ])
                    if bbr:
                        data_proc.connect([(bids_dg, pre_wf, [(('T1w_files', remove_container), 'coregwf.reg_bbr.reference')]),
                                          (pre_wf, datasink, [('coregwf.reg_bbr.out_matrix_file', 'preproc.' + pipelineID + '.@mat_file')]),
                                      ])
                    
                    
                        
                    l1 = self.first_level_flow(TR)
                    
                    data_proc.connect([(bids_dg, event_node, [(('bold_files', remove_container), 'file')]),
                                        (event_node, l1, [('event_file', 'modelspec.bids_event_file')]),
                                        (pre_wf, l1, [('art.outlier_files', 'modelspec.outlier_files')]),
                                        (pre_wf, l1, [('mc.par_file', 'modelspec.realignment_parameters')]),
                                        ])
                    
                    
                    if susan:
                        data_proc.connect(pre_wf, 'smooth.io_o.out_file', datasink, 'preproc.' + pipelineID + '.@smooth')
                        data_proc.connect(pre_wf, 'smooth.io_o.out_file', l1, 'modelspec.functional_runs')
                    else:
                        data_proc.connect(pre_wf, 'smooth.out_file', datasink, 'preproc.' + pipelineID + '.@smooth')
                        data_proc.connect(pre_wf, 'smooth.out_file', l1, 'modelspec.functional_runs')
                        
                    
                    norm = self.spatial_normalization_flow()
                    data_proc.connect([(pre_wf, norm, [('coregwf.bet_anat.out_file', 'warp.in_file')]),
                                       (pre_wf, norm, [('coregwf.bet_anat.out_file', 'prelim.in_file')]),
                                       (l1, norm, [('feat.feat_dir', 'selectfiles.base_directory')]),
                                       ])
                    
                    data_proc.run()
                    #TODO: SPLIT HALF STUFF
                    #TODO: FIRST LEVEL ANALYSIS
                    break
                else:
                    print("Error: Missing functional or structural images")
                        
        
        
    def genetic_algorithm(self, data):
        placeholder = True
        #TODO
        
    def ga_coder(self, decode=True):
        placeholder = True
    #INPUT OUTPUT NEXT
        
    def workflow_assembly(self, options, fwhm=4, iso_resample = 4, HP = 50, sig_loss_thresh = 20, bright_thresh = 2):
        
        
        #TR = 5
        #MAKE CONNECTIONS HERE
        A=3
        
    def spatial_normalization_flow(self):
        from nipype.interfaces.fsl import FNIRT, ApplyWarp, FLIRT
        from nipype import Node, MapNode, Workflow, SelectFiles
        import os
        from os.path import join as opj
        #COULD PROBABLY CHANGE THESE PARAMETERS - ASK ERIN, TONS FOR FNIRT
        
        ref_file = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
        prelim = Node(FLIRT(reference=ref_file, dof=12, output_type='NIFTI_GZ'), name='prelim')
        #WARP IS ACTING WEIRDLY, NOT OUTPUTTING WARP FILE
        warp = Node(FNIRT(ref_file=ref_file, field_file=True), name='warp')
        
        applywarp_t = MapNode(ApplyWarp(ref_file=ref_file), name='applywarp_t', iterfield=['in_file'])
        applywarp_z = MapNode(ApplyWarp(ref_file=ref_file), name='applywarp_z', iterfield=['in_file'])
        applywarp_c = MapNode(ApplyWarp(ref_file=ref_file), name='applywarp_c', iterfield=['in_file'])
        applywarp_v = MapNode(ApplyWarp(ref_file=ref_file), name='applywarp_v', iterfield=['in_file'])
        
        templates = {'t_stat': 'stats/tstat*.nii.gz',
                     'z_stat': 'stats/zstat*.nii.gz',
                     'cope': 'stats/cope*.nii.gz',
                     'varcope': 'stats/varcope*.nii.gz'}
        
        selectfiles = Node(SelectFiles(templates, sort_filelist=True), name='selectfiles')
        
        warpflow = Workflow('warpflow')
        warpflow.base_dir = opj(self.exp_dir, self.working_dir)
        
        warpflow.connect([(prelim, warp, [('out_matrix_file', 'affine_file')]),
                          (warp, applywarp_t, [('field_file', 'field_file')]),
                          (warp, applywarp_z, [('field_file', 'field_file')]),
                          (warp, applywarp_c, [('field_file', 'field_file')]),
                          (warp, applywarp_v, [('field_file', 'field_file')]),
                          (selectfiles, applywarp_t, [('t_stat', 'in_file')]),
                          (selectfiles, applywarp_z, [('z_stat', 'in_file')]),
                          (selectfiles, applywarp_c, [('cope', 'in_file')]),
                          (selectfiles, applywarp_v, [('varcope', 'in_file')]),
                          ])
        
        return warpflow
    
        
    def first_level_flow(self, TR, discard=4, HP=128, thresh=1000, serial_cor=True, base_switch=False, gamma=False, dict_opt='gammasigma', base_val=3):
        #MELODIC ICA DENOISING?
        #HIGH PASS FILTERING
        #FILM PREWHITENING (on/off)
        #MCFLIRT MOTION PARAMETERS AS CONFOUNDS EVs (explanatory variable) in model (either or with melodic ica denoising)
        #REMOVE MOTION RELATED OUTLIERS -> fsl_motion_outliers
        #PHYSIOLOGICAL NOISE MODEL -> probably not
        #MODEL SETUP FOR FILM
        #WHAT CONTRASTS DO WE WANT?
        from os.path import join as opj
        from nipype.interfaces.fsl import (FEATModel, FEATRegister, FEAT, FILMGLS, GLM, Level1Design)
        from nipype.algorithms.modelgen import SpecifyModel
        from nipype import Node, Workflow, Function
        
        #from os.path import join as opj
        #import json
        #from nipype.interfaces.spm import Level1Design, EstimateModel, EstimateContrast
        #from nipype.algorithms.modelgen import SpecifySPMModel
        #from nipype.interfaces.utility import Function, IdentityInterface
        #from nipype.interfaces.io import SelectFiles, DataSink
        #from nipype import Workflow, Node
        #import nipype.interfaces.matlab as matlab
        #from nipype.interfaces import spm
        #spm.SPMCommand.set_mlab_paths(matlab_cmd='/Applications/MATLAB_R2017b.app/bin/matlab')
        
        #matlab.MatlabCommand.set_default_matlab_cmd('/Applications/MATLAB_R2017b.app/bin/matlab')
        
        #Generates model
        #Outputs session_info
        #TO CONNECT: functional_runs, bids_event_file or event_files or subject_info, maybe outlier_files (art outlier files), maybe realignment parameters from MC
        modelspec = Node(SpecifyModel(input_units='secs', time_repetition=TR, high_pass_filter_cutoff=HP, 
                                      parameter_source='FSL'), name='modelspec')
# =============================================================================
#         modelspec = Node(SpecifySPMModel(concatenate_runs=False,
#                                  input_units='secs',
#                                  output_units='secs',
#                                  time_repetition=TR,
#                                  high_pass_filter_cutoff=128),
#                  name="modelspec")
#         
#         level1design = Node(Level1Design(bases={'hrf': {'derivs': [1, 0]}},
#                                  timing_units='secs',
#                                  interscan_interval=TR,
#                                  model_serial_correlations='FAST'),
#                     name="level1design")
#         
#         level1estimate = Node(EstimateModel(estimation_method={'Classical': 1}),
#                       name="level1estimate")
#         
#         level1conest = Node(EstimateContrast(), name="level1conest")
#         
#         l1analysis = Workflow(name='l1analysis')
#         l1analysis.base_dir = opj(self.exp_dir, self.working_dir)
# =============================================================================
        
        
        def generate_contrasts(session_info):
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
        
        def correct_task_timing(session_info, correction): #, time_correction=TR*discard):
            for i, info in enumerate(session_info):
                for j, task in enumerate(info['cond']):
                    for k, num in enumerate(task['onset']):
                        session_info[i]['cond'][j]['onset'][k] = num - correction
                        
            return session_info
        
        #def film_command_string(brightness_threshold, mask_size, autocorr_estimate_only=True, autocorr_noestimate=False, 
        #                        fit_armodel=False, multitaper_product=10, smooth_autocorr=True, 
        #                        threshold=1000, tukey_window=10, use_pava=False):
        # POTENTIAL TO ADD, PEOPLE GENERALLY DON'T CHANGE THESE IT SEEMS
            
        
        contrasts = Node(Function(input_names=['session_info'], output_names=['identities', 'contrasts'], 
                                  function=generate_contrasts), name='contrasts')
        
        correction = Node(Function(input_names=['session_info', 'correction'], output_names=['session_info'], 
                                  function=correct_task_timing), name='correction')
        
        correction.inputs.correction = TR * discard
        
        #l1analysis.connect([(modelspec, contrasts, [('session_info', 'session_info')]),
                            
        #                    (contrasts, level1conest, [('contrasts', 'contrasts')]),
        #            (modelspec, level1design, [(('session_info', correct_task_timing),
        #                                        'session_info')]),
        #            (level1design, level1estimate, [('spm_mat_file',
        #                                             'spm_mat_file')]),
        #            (level1estimate, level1conest, [('spm_mat_file',
        #                                             'spm_mat_file'),
        #                                            ('beta_images',
        #                                             'beta_images'),
        #                                            ('residual_image',
        #                                             'residual_image')]),
        #            ])
        
        #Generate EV files
        #basis_func options: dgamma
        #base_opt options: derivs(bool, or: dictionary gamma(derivs, gammasigma, gammadelay)
        #TO CONNECT: modelspec session_info output
        #OPTIONAL: contrasts, orthogonalization
        #OUTPUTS: ev_files, fsf_files
        if gamma:
            l1d = Node(Level1Design(bases={'gamma': {dict_opt: base_val}}, 
                                    interscan_interval=TR, 
                                    model_serial_correlations=serial_cor), name='l1d')
        else:
            l1d = Node(Level1Design(bases={'dgamma':{'derivs':base_switch}}, 
                                    interscan_interval=TR, 
                                    model_serial_correlations=serial_cor), name='l1d')
        
        #Generate design files
        #INPUTS: ev_files, fsf_files
        #OUTPUTS: design_file, con_file
        featmod = Node(FEATModel(output_type='NIFTI_GZ'), name='featmod')
        
        #Run FEAT
        feat = Node(FEAT(), name='feat')
        
        l1 = Workflow(name='l1')
        l1.base_dir = opj(self.exp_dir, self.working_dir)
        
        l1.connect([(modelspec, correction, [('session_info', 'session_info')]),
                    (correction, l1d, [('session_info', 'session_info')]),
                    (modelspec, contrasts, [('session_info', 'session_info')]),
                    (contrasts, l1d, [('contrasts', 'contrasts')]),
                    (l1d, featmod, [('fsf_files', 'fsf_file')]),
                    (l1d, featmod, [('ev_files', 'ev_files')]),
                    (l1d, feat, [('fsf_files', 'fsf_file')]),
                    #(featmod, feat, [('design_file', 'fsf_file')]),
                    ])
        
        
        #WILL LIKELY NEED TO GENERATE CONTRASTS
        
        #ABSURD LIST OF OPTIONAL INPUTS
        #IN: input data file
        #OPT: design_file, brightness_thresh & mask_size (from SUSAN)
        #OUT: corrections, dof_file, param_estimates, sigmasquareds, thresholdac
        #film = Node(FILMGLS(threshold=thresh, autocorr_noestimate=not serial_cor, output_type='NIFTI_GZ'), name='film')
        
        return l1
        
        
        
        
        
        
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
            
            
    
    #might want to add non-linear registration as an option FNIRT
    def coreg_flow(self, options, bet_frac=0.5, wm_thresh=0.5, dof_f=6, bbr=True, bbr_type='signed', interp='spline', iso=4, cost='mutualinfo', bins=640): #iso -> iso_resample
        #NOTE: this code is adapted from the nipype tutorial on preprocessing https://miykael.github.io/nipype_tutorial/notebooks/example_preprocessing.html 
        from os.path import join as opj
        import os
        from nipype.interfaces.fsl import (BET, FAST, FLIRT, Threshold)
        from nipype import Workflow, Node, MapNode
        
        #Brain extraction (maybe frac is a parameter that can be adjusted, could alter -g vertical gradient intensity threshold)
        bet_anat = Node(BET(frac=bet_frac, robust=True, output_type='NIFTI_GZ'), name="bet_anat")#, iterfield='in_file')
        #FAST image segmentation - necessary for bbr registration (recommended) -> get white matter segmentation
        segment = Node(FAST(output_type='NIFTI_GZ'), name='segment')
        #Get white matter segmentation output
        def get_wm(files):
            return files[-1]
        
        #threshold white matter probability mask to get binary image -> can adjust thresh
        threshold = Node(Threshold(thresh=wm_thresh, args='-bin', output_type='NIFTI_GZ'), name="threshold")
        #flirt - register functional to anatomical images -> vanilla (no bbr)
        #cost=cost
        reg_pre = Node(FLIRT(dof=dof_f, output_type='NIFTI_GZ'), name='reg_pre')
        #flirt - bbr registration (there are 3 different types of bbr)
        #schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch'), Doesn't exist
        #bbrtype=bbr_type
        reg_bbr = Node(FLIRT(dof=dof_f, cost='bbr', schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch'), output_type='NIFTI_GZ'), name='reg_bbr')
        #apply registration warp to functional images
        applywarp = Node(FLIRT(interp=interp, apply_isoxfm=iso, output_type='NIFTI_GZ'), name='applywarp')
        #apply registration warp to mean file (seems to be the same as above just gunzipped - no smoothing)
        applywarp_mean = Node(FLIRT(interp=interp, apply_isoxfm=iso, output_type='NIFTI_GZ'), name='applywarp_mean')
        
        #coregistration workflow
        coregwf = Workflow(name='coregwf')
        coregwf.base_dir = opj(self.exp_dir, self.working_dir)
        
        #node connection
        coregwf.connect([(bet_anat, reg_pre, [('out_file', 'reference')]),
                         (bet_anat, applywarp, [('out_file', 'reference')]),
                         (bet_anat, applywarp_mean, [('out_file', 'reference')]),
                         ])
        if bbr:
            coregwf.connect([(bet_anat, segment, [('out_file', 'in_files')]),
                             (segment, threshold, [(('partial_volume_files', get_wm), 'in_file')]),
                             (threshold, reg_bbr, [('out_file', 'wm_seg')]),
                             (reg_pre, reg_bbr, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, applywarp, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_bbr, applywarp_mean, [('out_matrix_file', 'in_matrix_file')]),
                             ])
        else:
            coregwf.connect([(reg_pre, applywarp, [('out_matrix_file', 'in_matrix_file')]),
                             (reg_pre, applywarp_mean, [('out_matrix_file', 'in_matrix_file')]),
                             ])
            
        return coregwf
            
        
A = fmri_multiverse()
        
        
        
        
        
        