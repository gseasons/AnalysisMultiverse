#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:50:42 2021

@author: grahamseasons
"""

#have so that multiple values inputted as tuples
#resting will be a dictionary with everything seed related inside
#multiple seeds, separate statistical maps

def info(mask, task, TR, event_file, unsmoothed, smoothed, outliers, brain, brainmask, segmentations, invwarp, network):
    from nipype import Node
    from versatile import SpecifyModelVersatile
    from nipype.interfaces.fsl import ImageMeants, ExtractROI, ImageMaths, ApplyWarp, Threshold, BinaryMaths, WarpPoints
    from nipype.interfaces.fsl.maths import MathsCommand
    from functions import parse_xml
    import re, os
    import numpy as np
    from updated.l1_analysis.functions import data_driven, warp, invert
    
# =============================================================================
#     warppostfeat=True
#     concatenate=True
#     #rest = {'type': 'ROI', 'coords': [(50, 25, 60), (75, 30, 60), (100, 85, 62)], 'radius': 8.0}
#     rest = {'type': 'data', 'seedinfo': [('cing post', 17)], 'atlas': 'harvard', 'k': 'faces', 'kcc': 0.7, 'lp': 0.01, 'hp': 0.1}
# =============================================================================
    
    model = Node(SpecifyModelVersatile(input_units='secs', parameter_source='FSL'), name='model')
    model.inputs.outlier_files = outliers
    model.inputs.time_repetition = TR
    model.inputs.high_pass_filter_cutoff = vars().get('HP', 128)
    model.inputs.functional_runs = smoothed
    
    if event_file:
        model.inputs.bids_event_file = event_file
        session_info = model.run().outputs.session_info
        
        return session_info
    elif 'rest' in task and vars().get('rest'):
        #NOTE: COORDINATES WILL BE IN STANDARD SPACE, AS WILL ATLAS LOCATIONS -> TRANSFER TO STANDARD TO GET TIME COURSES
        if vars()['rest']['type'] == 'ROI':
            x, y, z = vars()['rest']['coords'][network]
            
            if vars()['warppostfeat']:
                roipoints = os.getcwd() + '/roipoints.txt'
                np.savetxt(roipoints, np.array([x,y,z]))
                if vars()['concatenate']:
                    img_file = brainmask
                    mask_ = brainmask
                else:
                    img_file = brain
                    mask_ = brain
                warppoints = WarpPoints(in_coords=roipoints, dest_file=img_file, src_file=mask, warp_file=invwarp).run().outputs.out_file
                points = np.loadtxt(warppoints)
                x = points[0]
                y = points[1]
                z = points[2]
            else:
                mask_ = mask
            
            radius = vars()['rest']['radius']
            createseed = Node(MathsCommand(in_file=mask_), name='createseed')
            createseed.inputs.args = '-mul 0 -add 1 -roi {x} 1 {y} 1 {z} 1 0 1'.format(x=x, y=y, z=z)
            seed = createseed.run().outputs.out_file
            
            makesphere = Node(MathsCommand(in_file=seed), name='makesphere')
            makesphere.inputs.args = '-kernel sphere {radius} -fmean'.format(radius=radius)
            
            thrseed = ImageMaths(op_string='-bin')
            thrseed.inputs.in_file = makesphere.run().outputs.out_file
            thrfile = thrseed.run().outputs.out_file
            suffix = '_x{x}_y{y}_z{z}_r{r}'.format(x=x, y=y, z=z, r=radius)
        elif vars()['rest']['type'] == 'atlas' or vars()['rest']['type'] == 'data':
            #nipype.algorithms.misc PickAtlas(although i've kind of implemented this, might be nice to add multiple ROIs -> get user input for options and index them in GA, pass list of list)
            #ADD WHEN SETTING UP USER INPUT GUI -> atlas will be link to atlas to be used by PickAtlas
            seed, thr = vars()['rest']['seedinfo'][network]
            atlas = vars()['rest']['atlas']
            file, index, name = parse_xml(atlas, seed, mask)
            getseed = Node(ExtractROI(in_file=file, t_min=int(index), t_size=1), name='getseed')
            roi = getseed.run().outputs.roi_file
            
            if vars()['warppostfeat']:
                if vars()['concatenate']:
                    ref_file = brainmask
                else:
                    ref_file = brain
                roi = ApplyWarp(in_file=roi, ref_file=ref_file, field_file=invwarp).run().outputs.out_file
            
            thrseed = ImageMaths(op_string='-thr {thr} -bin'.format(thr=thr))
            thrseed.inputs.in_file = roi
            thrfile = thrseed.run().outputs.out_file
            suffix = ''.join([word[0] for word in name.split()])
            
            if vars()['rest']['type'] == 'data':
                k = vars()['rest']['k']
                kcc = vars()['rest']['kcc']
                lp = vars()['rest'].get('lp', 0.01)
                hp = vars()['rest'].get('hp', 0.1)
                
                reho = data_driven(ref_file, unsmoothed, k, kcc, TR, lp, hp)
                mul = BinaryMaths(in_file=reho, operand_file=thrfile, operation='mul')
                thrfile = mul.run().outputs.out_file
                suffix += '_reho'
        else:
            print('ERROR')
            
        ev_name = re.search('task-([a-zA-Z]+)_', smoothed).group(1) + suffix
                    
        mean_ts = Node(ImageMeants(in_file=smoothed, out_file=ev_name), name='mean_ts')
                    
        mean_ts.inputs.mask = thrfile
        time_series = [mean_ts.run().outputs.out_file]
                
        model.inputs.event_files = time_series
        session_info = model.run().outputs.session_info
        
        return session_info, thrfile
    else:
        print('ERROR')
    

# =============================================================================
# TR=3
# brain = '/Volumes/NewVolume/sup_pre_l1/working_dir/fmri/preprocess/bet/atropos_wf/_subject_002S6009/apply_mask/mapflow/_apply_mask0/sub-002S6009_T1w_corrected_xform_masked.nii.gz'
# brainmask = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/decision/sub-002S6009_T1w_corrected_xform_masked_thresh_flirt.nii.gz'
# event_file = ''
# invwarp = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/invwarp/sub-002S6009_T1w_corrected_xform_masked_thresh_flirt_concatwarp.nii.gz'
# mask = '/usr/local/fsl/data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz'
# network = 0
# outliers = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/art/art.sub-002S6009_task-rest_bold_roi_mcf_outliers.txt'
# segmentations = ['/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_0_flirt.nii.gz', '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_1_flirt.nii.gz', '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_2_flirt.nii.gz']
# smoothed = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/_i_0/_i_0/Fsmooth/smooth/smooth_su/sub-002S6009_task-rest_bold_roi_mcf_regressed_maths_smooth.nii.gz'
# task = 'rest'
# unsmoothed = '/Volumes/NewVolume/sup_pre_l1_v2/working_dir/fmri/preprocess/_subject_002S6009/Fregress/sub-002S6009_task-rest_bold_roi_mcf_regressed_maths.nii.gz'
# 
# info(mask, task, TR, event_file, unsmoothed, smoothed, outliers, brain, brainmask, segmentations, invwarp, network)
# =============================================================================

