#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:48:57 2021

@author: grahamseasons
"""
from niworkflows.anat.ants import init_brain_extraction_wf
from bids.layout import BIDSLayout
from nipype import IdentityInterface, Node, DataSink
from nipype import config, Workflow, Node, IdentityInterface
#config.set("execution", "remove_node_directories", "true")
# =============================================================================
# global _10, _67, _3
# _10 = [55, 67]
# _67 = [56, 1,2,3,4]
# _3 = 4502
# 
# 
# def _graph_expansion(lst, searched):
#         lst_out = lst.copy()
#         for idx in lst:
#             if idx in searched:
#                 continue
#             lst_out += vars().get('_'+str(idx), [])
#             searched.append(idx)
#         
#         return len(searched) != len(lst_out), lst_out, searched
#     
# lst = [10]
# searched = []
# egg = True
# while egg:
#     egg, lst, searched = _graph_expansion(lst, searched)
# =============================================================================

def info(mask, task, TR, event_file, unsmoothed, smoothed, brain, brainmask, outliers, segmentations, invwarp, network, rest, HP, warppostfeat, concatenate):
    from nipype import Node
    from versatile import SpecifyModelVersatile
    from nipype.interfaces.fsl import ImageMeants, ExtractROI, ImageMaths, BinaryMaths, WarpPoints
    from nipype.interfaces.fsl.maths import MathsCommand
    from l1_analysis.functions import parse_xml
    import re, os
    import numpy as np
    from l1_analysis.functions import data_driven, warp

    model = Node(SpecifyModelVersatile(input_units='secs', parameter_source='FSL'), name='model')
    model.inputs.time_repetition = TR
    model.inputs.high_pass_filter_cutoff = vars().get('HP', 128)
    model.inputs.functional_runs = smoothed
    model.inputs.outlier_files = outliers

    if event_file:
        model.inputs.bids_event_file = event_file
        session_info = model.run().outputs.session_info
        session_info[0]['scans'] = smoothed
        return session_info
    elif 'rest' in task and vars().get('rest'):
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
                points = np.loadtxt(warppoints).astype(int)
                x = points[0]
                y = points[1]
                z = points[2]
            else:
                mask_ = mask

            radius = vars()['rest']['radius']
            createseed = Node(MathsCommand(in_file=mask_, output_datatype='float'), name='createseed')
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
            if isinstance(vars()['rest']['seedinfo'][network], (tuple, list)):
                seed, thr = vars()['rest']['seedinfo'][network]
                atlas = vars()['rest']['atlas']
                file, index, name = parse_xml(atlas, seed, mask)
                suffix = ''.join([word[0] for word in name.split()])
                getseed = Node(ExtractROI(in_file=file, t_min=int(index), t_size=1), name='getseed')
                roi = getseed.run().outputs.roi_file
            else:
                roi = vars()['rest']['seedinfo'][network]
                suffix = re.search('.*([A-Za-z0-9_-]+).nii.*').group(1)
                thr = 0

            if vars()['warppostfeat']:
                if vars()['concatenate']:
                    ref_file = brainmask
                else:
                    ref_file = brain
                roi = warp(roi, ref_file, invwarp)
            else:
                ref_file = brainmask

            thrseed = ImageMaths(op_string='-thr {thr} -bin'.format(thr=thr))
            thrseed.inputs.in_file = roi
            thrfile = thrseed.run().outputs.out_file


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
            raise NotImplementedError("Invalid seeding method of {method} used, which is not implemented. Please use 'atlas', 'ROI', or 'data'.".format(method=vars()['rest']['type']))

        ev_name = re.search('task-([a-zA-Z]+)_', smoothed).group(1) + suffix

        mean_ts = Node(ImageMeants(in_file=smoothed, out_file=ev_name), name='mean_ts')

        mean_ts.inputs.mask = thrfile
        time_series = [mean_ts.run().outputs.out_file]

        model.inputs.event_files = time_series
        session_info = model.run().outputs.session_info

        session_info[0]['scans'] = smoothed

        return session_info, thrfile
    else:
        raise ValueError("Unhandled task of {task}. If resting state analysis ensure 'rest' is in the task name, otherwise ensure there is a valid event_file".format(task=task))

mask = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz'
task = 'rest'
TR = 3
event_file=''
unsmoothed = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_0/_subject_002S6009/_i_5/_i_5/_i_5/Fregress/sub-002S6009_task-rest_bold_roi_mcf_st_warp_reho_regressed_maths.nii.gz'
smoothed = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_0/_subject_002S6009/_i_5/_i_5/_i_5/Fsmooth/smooth/smooth_su/sub-002S6009_task-rest_bold_roi_mcf_st_warp_regressed_maths_smooth.nii.gz'
brain = '/Volumes/NewVolumes/super_agers/brain_extracted/_subject_002S6009/sub-002S6009_T1w_corrected_xform_masked.nii.gz'
brainmask = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_0/_subject_002S6009/_i_5/_i_5/_i_5/fillmask/uni_xform_masked_fillh.nii.gz'
outliers = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_subject_002S6009/_i_0/_i_5/_i_5/_i_5/art/art.sub-002S6009_task-rest_bold_roi_mcf_st_outliers.txt'
segmentations = ['/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_0/_subject_002S6009/_i_5/_i_5/_i_5/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_0_warp.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_0/_subject_002S6009/_i_5/_i_5/_i_5/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_1_warp.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_0/_subject_002S6009/_i_5/_i_5/_i_5/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_2_warp.nii.gz']
invwarp = ''
network = 0
rest = {'type': 'data', 'seedinfo': [('Cingulate Gyrus, posterior division', 1)], 'atlas': 'harvard', 'k': 'faces', 'kcc': 0.6, 'lp': 0.01, 'hp': 0.1}
HP = 128
warppostfeat = False
concatenate = True
#info(mask, task, TR, event_file, unsmoothed, smoothed, brain, brainmask, outliers, segmentations, invwarp, network, rest, HP, warppostfeat, concatenate)


def mni(mniMask, brain, brainmask, warp, warped, segmentations, out_mat, start_img, warplater, warped_interp, concatenate):
    from nipype.interfaces.fsl import ApplyWarp, ConvertWarp, Threshold, ConvertXFM, FLIRT
    import re
    from nipype.interfaces.base import Undefined
    if not vars().get('warplater', True):

        for i, file in enumerate(segmentations):
            segmentations[i] = ApplyWarp(in_file=file, ref_file=mniMask, field_file=warp, interp='nn').run().outputs.out_file

        if vars().get('concatenate', True):
            warp = ConvertWarp(reference=mniMask, premat=out_mat, warp1=warp).run().outputs.out_file

        warped = ApplyWarp(in_file=warped, ref_file=mniMask, field_file=warp)

        for param in ['warped_interp']:
                search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
                if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
                else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

        warped = warped.run().outputs.out_file

        start_warped = ApplyWarp(in_file=start_img, ref_file=mniMask, field_file=warp)

        start_img = start_warped.run().outputs.out_file
        brainmask = ApplyWarp(in_file=brainmask, ref_file=mniMask, field_file=warp, interp='nn').run().outputs.out_file
    else:
        if not vars().get('concatenate', False):
            brainmask = Threshold(in_file=brain, thresh=0, args='-bin').run().outputs.out_file
        else:
            warp = ConvertWarp(reference=mniMask, premat=out_mat, warp1=warp).run().outputs.out_file
            segform = ConvertXFM(in_file=out_mat, invert_xfm=True).run().outputs.out_file
            for i, file in enumerate(segmentations):
                segmentations[i] = FLIRT(in_file=file, reference=warped, interp='nearestneighbour', apply_xfm=True, in_matrix_file=segform).run().outputs.out_file

    return warped, brainmask, segmentations, warp, start_img

def regress(unsmoothed, mc_par, segmentations, mask, rest, CSF, WM, GLOBAL, glm_dat_norm, glm_demean, realignregress, glm_des_norm):
    from nipype import Node
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.fsl import ImageMeants, Threshold, FLIRT, GLM, ImageStats, ImageMaths, ExtractROI
    from nipype.interfaces.fsl.maths import MeanImage
    from nipype.interfaces.base import CommandLine
    import numpy as np
    import nibabel as nib
    import os, re
    base_dir = os.getcwd()
    
    if rest:
        CSF = vars().get('CSF', True)
        WM = vars().get('WM', True)
        GLOBAL = vars().get('GLOBAL', True)
    else:
        CSF = vars().get('CSF', False)
        WM = vars().get('WM', False)
        GLOBAL = vars().get('GLOBAL', False)

    params = np.loadtxt(mc_par)
    reho = np.loadtxt(mc_par)

    if not vars().get('realignregress', True):
        params = np.zeros((params.shape[0], 1))

    resample = False
    suffix = ''

    reference = Node(MeanImage(in_file=unsmoothed, dimension='T'), name='reference', base_dir=base_dir).run().outputs.out_file
    if nib.load(unsmoothed).shape[0:-1] != nib.load(mask).shape:
        resample = True

    csfmask = segmentations[0]
    csfmask = ImageMaths(in_file=csfmask, in_file2=mask, op_string='-mul').run().outputs.out_file
    if resample:
        csfmask = Node(FLIRT(in_file=csfmask, reference=reference, apply_xfm=True, uses_qform=True), name='csfmask', base_dir=base_dir).run().outputs.out_file
    meancsf = Node(ImageMeants(in_file=unsmoothed, mask=csfmask), name='meancsf', base_dir=base_dir)
    csf = np.loadtxt(meancsf.run().outputs.out_file).reshape(-1, 1)
    reho = np.hstack((reho, csf))

    if CSF:
        params = np.hstack((params, csf))
        suffix += 'CSF_'

    wmmask = segmentations[2]
    wmmask = ImageMaths(in_file=wmmask, in_file2=mask, op_string='-mul').run().outputs.out_file
    if resample:
        wmmask = Node(FLIRT(in_file=wmmask, reference=reference, apply_xfm=True, uses_qform=True), name='wmmask', base_dir=base_dir).run().outputs.out_file
    meanwm = Node(ImageMeants(in_file=unsmoothed, mask=wmmask), name='meancsf', base_dir=base_dir)
    wm = np.loadtxt(meanwm.run().outputs.out_file).reshape(-1, 1)
    reho = np.hstack((reho, wm))

    if WM:
        params = np.hstack((params, wm))
        suffix += 'WM_'

    if GLOBAL:
        meanglob = Node(ImageMeants(in_file=unsmoothed, mask=mask), name='meanglob', base_dir=base_dir)
        glob = np.loadtxt(meanglob.run().outputs.out_file).reshape(-1, 1)
        params = np.hstack((params, glob))
        reho = np.hstack((reho, glob))
        suffix += 'GLOBAL'

    name_ = os.getcwd() + '/' + suffix

    if not vars().get('realignregress', True):
        params = params[:,1:]

    np.savetxt(name_ + '.txt', params)
    cmd = ('Text2Vest {name_}.txt {name_}.mat').format(name_=name_)
    cl = CommandLine(cmd)
    cl.run().runtime

    name_reho = os.getcwd() + '/' + suffix + '_reho'
    np.savetxt(name_ + '.txt', reho)
    cmd = ('Text2Vest {name_}.txt {name_reho}.mat').format(name_=name_, name_reho=name_reho)
    cl = CommandLine(cmd)
    cl.run().runtime

    forreho = unsmoothed
    if np.any(params):
        out_name = re.search('/([A-Za-z0-9_-]+).nii', unsmoothed).group(1) + '_regressed.nii.gz'
        glm = GLM(design=name_ + '.mat', in_file=unsmoothed, out_res_name=out_name)

        for param in ['glm_dat_norm', 'glm_demean', 'glm_des_norm']:
                search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
                if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
                else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

        regressed = glm.run().outputs.out_res
        add = abs(ImageStats(in_file=regressed, op_string='-R').run().outputs.out_stat[0])
        unsmoothed = ImageMaths(in_file=regressed, args='-add '+str(add)).run().outputs.out_file

    if np.any(reho) and not np.array_equal(reho, params):
        out_name = re.search('/([A-Za-z0-9_-]+).nii', forreho).group(1) + '_reho_regressed.nii.gz'
        glm = GLM(design=name_reho + '.mat', in_file=forreho, out_res_name=out_name)

        for param in ['glm_dat_norm', 'glm_demean', 'glm_des_norm']:
                search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
                if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
                else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

        regressed = glm.run().outputs.out_res
        add = abs(ImageStats(in_file=regressed, op_string='-R').run().outputs.out_stat[0])
        forreho = ImageMaths(in_file=regressed, args='-add '+str(add)).run().outputs.out_file

    return unsmoothed, forreho

def smooth(warped, mask, susan, fwhm):
    #WORKFLOW ADAPTED FROM: https://nipype.readthedocs.io/en/latest/users/examples/fmri_fsl.html
    from nipype import Workflow, Node, IdentityInterface
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.fsl import (ImageMaths, ImageStats, SUSAN, IsotropicSmooth)
    import os, glob, re
    from preprocessing.functions import get_bright_thresh, getthreshop

    smooth = Workflow('smooth')
    smooth.base_dir = os.getcwd()

    susan = vars().get('susan', True)
    fwhm = vars().get('fwhm', 6)

    inputnode = Node(IdentityInterface(fields=['in_file']), name='inputnode')
    inputnode.inputs.in_file = warped

    outnode = Node(IdentityInterface(fields=['smoothed']), name='outnode')

    if susan:
        maskfunc = Node(ImageMaths(suffix='_bet', op_string='-mas', in_file2=mask), name='maskfunc')
        getthresh = Node(ImageStats(op_string='-p 2 -p 98'), name='getthreshold')
        threshold = Node(ImageMaths(out_data_type='char', suffix='_thresh'), name='threshold')

        medianval = Node(ImageStats(op_string='-k %s -p 50'), name='medianval')

        smooth_su = Node(SUSAN(fwhm=fwhm), name='smooth_su')

        smooth.connect([(inputnode, maskfunc, [('in_file', 'in_file')]),
                        (maskfunc, getthresh, [('out_file', 'in_file')]),
                        (maskfunc, threshold, [('out_file', 'in_file')]),
                        (getthresh, threshold, [(('out_stat', getthreshop), 'op_string')]),
                        (inputnode, medianval, [('in_file', 'in_file')]),
                        (threshold, medianval, [('out_file', 'mask_file')]),
                        (medianval, smooth_su, [(('out_stat', get_bright_thresh), 'brightness_threshold')]),
                        (inputnode, smooth_su, [('in_file', 'in_file')]),
                        (smooth_su, outnode, [('smoothed_file', 'smoothed')]),
                        ])

        for param in []:
                search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
                if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
                else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

        smooth.run()
        smoothed = glob.glob(os.getcwd() + '/smooth/smooth_su/*.nii*')[0]
    else:
        smooth_iso = Node(IsotropicSmooth(fwhm=fwhm, output_type='NIFTI'), name='smooth_iso')
        smooth.connect(inputnode, 'in_file', smooth_iso, 'in_file')
        smooth.connect(smooth_iso, 'out_file', outnode, 'smoothed')

        for param in []:
                search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)
                if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])
                else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)

        smooth.run()
        smoothed = glob.glob(os.getcwd() + '/smooth/smooth_iso/*.nii*')[0]

    return smoothed, glob.glob(os.getcwd() + '/smooth/**', recursive=True)

fwhm = 12
susan = True
warped = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/Fregress/sub-002S6009_task-rest_bold_roi_mcf_st_warp_regressed_maths.nii.gz'
mask = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/fillmask/uni_xform_masked_fillh.nii.gz'
#smooth(warped, mask, susan, fwhm)

segmentations = ['/Users/grahamseasons/fMRI/analysis_multiverse/hide/testing/code/l1_analysis/sub-002S6009_T1w_corrected_xform_masked_seg_2_flirt_warp.nii.gz',
 '/Users/grahamseasons/fMRI/analysis_multiverse/hide/testing/code/l1_analysis/sub-002S6009_T1w_corrected_xform_masked_seg_0_flirt_warp.nii.gz',
 '/Users/grahamseasons/fMRI/analysis_multiverse/hide/testing/code/l1_analysis/sub-002S6009_T1w_corrected_xform_masked_seg_1_flirt_warp.nii.gz']#['/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_subject_002S6009/fast/sub-002S6009_T1w_corrected_xform_masked_seg_0.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_subject_002S6009/fast/sub-002S6009_T1w_corrected_xform_masked_seg_1.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_subject_002S6009/fast/sub-002S6009_T1w_corrected_xform_masked_seg_2.nii.gz']#['/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_0_warp.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_1_warp.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/Fmni/sub-002S6009_T1w_corrected_xform_masked_seg_2_warp.nii.gz']
unsmoothed = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/Fmni/sub-002S6009_task-rest_bold_roi_mcf_st_warp.nii.gz'
mc_par = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_subject_002S6009/_i_2/_i_2/mcflirt/sub-002S6009_task-rest_bold_roi_mcf.nii.gz.par'
mask = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/fillmask/uni_xform_masked_fillh.nii.gz' #'/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/Fmni/sub-002S6009_T1w_corrected_xform_masked_thresh_flirt_warp.nii.gz'
rest = True
realignregress = False
CSF = True
WM = True
GLOBAL = True
glm_dat_norm = True
glm_demean = True
glm_des_norm = True
#regress(unsmoothed, mc_par, segmentations, mask, rest, CSF, WM, GLOBAL, glm_dat_norm, glm_demean, realignregress, glm_des_norm)

mniMask = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz'
brain = '/Volumes/NewVolume/super_agers/brain_extracted/_subject_002S6009/sub-002S6009_T1w_corrected_xform_masked.nii.gz'
brainmask = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_subject_002S6009/_i_2/_i_2/_i_2/decision/sub-002S6009_T1w_corrected_xform_masked_thresh_flirt.nii.gz'
warp = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_i_2/_subject_002S6009/_i_2/_i_2/_i_2/warp/sub-002S6009_T1w_corrected_xform_masked_field.nii.gz'
warped = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_subject_002S6009/_i_2/_i_2/slicetimer/sub-002S6009_task-rest_bold_roi_mcf_st.nii.gz'
segmentations = ['/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_subject_002S6009/fast/sub-002S6009_T1w_corrected_xform_masked_seg_0.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_subject_002S6009/fast/sub-002S6009_T1w_corrected_xform_masked_seg_1.nii.gz', '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_subject_002S6009/fast/sub-002S6009_T1w_corrected_xform_masked_seg_2.nii.gz']
out_mat = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_subject_002S6009/_i_2/_i_2/_i_2/_i_2/Fregistration/reg/regbbr/sub-002S6009_task-rest_bold_roi_mcf_roi_flirt.mat'
start_img = '/Volumes/NewVolume/AANGST/working_dir/fmri/preprocess/_i_0/_subject_002S6009/_i_2/_i_2/_i_2/decision/sub-002S6009_task-rest_bold_roi_mcf_roi.nii.gz'
warplater = False
warped_interp = 'trilinear'
concatenate = True
#mni(mniMask, brain, brainmask, warp, warped, segmentations, out_mat, start_img, warplater, warped_interp, concatenate)

import ipyparallel, os
data_dir = '/Volumes/NewVolume/super_agers'
out_dir = '/Volumes/NewVolume/brain_extracted'
layout = BIDSLayout(data_dir)
subjects = layout.get_subjects()
data = Node(DataSink(), name='data')
config.set("execution", "remove_node_directories", "true")
cluster = ipyparallel.Cluster(n=2)
cluster.start_cluster_sync()
for subject in subjects:
    wf = init_brain_extraction_wf()
    wf.base_dir = os.getcwd()
    wf.inputs.inputnode.in_files = '/Volumes/NewVolume/super_agers/sub-{SUB}/anat/sub-{SUB}_T1w.nii.gz'.format(SUB=subject)
    #wf.connect(wf.get_node('outputnode'), 'out_file', data, out_dir)
    wf.run(plugin='IPython', plugin_args={'profile': 'default', 'cluster_id': cluster.cluster_id})

from nipype.algorithms.modelgen import orth
import numpy as np
from numpy.linalg import lstsq
from scipy.linalg import orth
import os
from nipype import Node, Workflow, IdentityInterface, Function, JoinNode
#from nipype.interaces import IdentityInterface


# =============================================================================
werk = Workflow('werk')
werk.base_dir = os.getcwd()
ident = Node(IdentityInterface(fields=['iter']), name='ident')
ident.iterables = [('iter', [0,4])]
# =============================================================================

def printer(inp):
    print(inp)
    
def ret(inp):
    A = ['A', 'B', 'C', 'D', 'E']
    return A[inp]

def printer2(inp, inp2):
    print(inp)
    print(inp2)

test = Node(IdentityInterface(fields=['func_str']), name='test')
test.inputs.func_str = 'def printer(inp):\n print(inp)'

def greeeg(greg):
    return greg

dummy = Node(Function(input_names='greg', output_names='greg', function=greeeg), name='dummy')

# =============================================================================
# =============================================================================
# ident2 = Node(IdentityInterface(fields=['dum', 'iter']), name='ident2')
# ident2.itersource = ('ident', 'iter')
# ident2.iterables = [('dum', {0:[1,2], 4:[4]})]
# =============================================================================
# =============================================================================
####MAKE MINIMUM EXAMPLE TO REPRODUCE NODES BEING REMOVED TOO EARLY -> DELETE ON LAST TIME USED NOT FIRST
f = Node(Function(input_names='inp', function=printer), name='f')
f2 = Node(Function(input_names='inp', function=printer), name='f2')
f3 = Node(Function(input_names='inp', output_names='inp', function=ret), name='f3')
f4 = Node(Function(input_names='inp', function=printer2), name='f4')
f.iterables = [('inp', {0: {1: [2,3], 4: [5,6]}})]
join = Node(IdentityInterface(fields=['dum', 'egg']), name='join')#, joinsource='ident2', joinfield=['dum', 'egg'])
join2 = JoinNode(IdentityInterface(fields=['dum', 'egg']), name='join2', joinsource='ident', joinfield=['dum', 'egg'])
# =============================================================================
werk.connect([#(test, ident, [('func_str', 'function_str')]),
              #(test, ident2, [('func_str', 'function_str')]),
              (ident, f2, [('iter', 'inp')]),
              #(ident, ident2, [('iter', 'iter')]),
              # (ident2, f, [('dum', 'inp')]),
              # (ident2, join, [('dum', 'dum')]),
               (ident, f3, [('iter', 'inp')]),
               (f3, dummy, [('inp', 'greg')]),
               (dummy, join, [('greg', 'egg')]),
               (join, join2, [('dum', 'dum'),
                              ('egg', 'egg')]),
               (join2, f4, [('dum', 'inp'),
                            ('egg', 'inp2')]),
               ])
# =============================================================================

werk.run(plugin='MultiProc')#plugin='MultiProc', plugin_args={'n_procs': 1, 'memory_gb': 0.5})
import sys
sys.exit()

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
        
copes = [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]
varcopes = [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']]
groups = [[['002S6030', '002S6053', '003S4644', '003S6014', '009S6212'], ['002S6009', '002S6103', '003S6067', '003S6256', '006S0731']], [['002S6009', '002S6030', '003S6014', '003S6067', '003S6256'], ['002S6053', '002S6103', '003S4644', '006S0731', '009S6212']], [['002S6030', '002S6053', '003S4644', '003S6067', '003S6256'], ['002S6009', '002S6103', '003S6014', '006S0731', '009S6212']], [['002S6009', '002S6030', '002S6053', '006S0731', '009S6212'], ['002S6103', '003S4644', '003S6014', '003S6067', '003S6256']], [['002S6009', '002S6053', '002S6103', '003S4644', '003S6256'], ['002S6030', '003S6014', '003S6067', '006S0731', '009S6212']], [['002S6009', '002S6030', '002S6103', '003S4644', '006S0731'], ['002S6053', '003S6014', '003S6067', '003S6256', '009S6212']], [['002S6009', '002S6053', '003S4644', '006S0731', '009S6212'], ['002S6030', '002S6103', '003S6014', '003S6067', '003S6256']], [['002S6009', '002S6103', '003S4644', '003S6014', '006S0731'], ['002S6030', '002S6053', '003S6067', '003S6256', '009S6212']], [['002S6103', '003S4644', '003S6014', '003S6256', '009S6212'], ['002S6009', '002S6030', '002S6053', '003S6067', '006S0731']], [['002S6053', '003S4644', '003S6014', '003S6256', '009S6212'], ['002S6009', '002S6030', '002S6103', '003S6067', '006S0731']]]
#construction(groups, copes, varcopes)

def find_orth(O, vec):
    A = np.hstack((O, vec))
    b = np.zeros(O.shape[1] + 1)
    b[-1] = 5
    return A, lstsq(A.T, b, rcond=None)[0]

A = {'normal': [1, 1, 1, 0, 0, 0, 0, 1, 0, 1], 'super': [0, 0, 0, 1, 1, 1, 1, 0, 1, 0], 'age': [67.6, 65.2, 65.7, 70.0, 72.8, 67.1, 63.1, 66.4, 78.9, 71.4], 'sex': [1, 0, 1, 0, 0, 1, 0, 0, 1, 0], 'education': [16, 18, 16, 16, 14, 18, 18, 16, 18, 18]}
B = np.array(list(A.values()))
test = np.array([[1,1,1,0,0,0,0,1,0,1], [0,0,0,1,1,1,1,0,1,0]]).T
#A, B = find_orth(test, np.array([1,1,1,0,0,1,0,1,0,1], ndmin=2).T)
#np.array()

A = np.array([[1,0,1,0,0,0,1], [0,1,0,1,1,1,0]])
#orth([1, 2, 3], [0, 1, 2])

def buffer(copes, varcopes, run_mode=[0,0], i=0, rng=2):
            if rng > 1:
                new_c = []
                new_v = []
                for k in range(len(copes)):
                    new_c.append(copes[k][i][0])
                    new_v.append(varcopes[k][i][0])
                    
                return new_c, new_v, run_mode[i]
            else:
                return copes[0], varcopes[0], run_mode
            
copes = [[[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]], [[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_0/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_1/applywarp_c/mapflow/_applywarp_c0/cope1_warp.nii.gz']]]]
varcopes = [[[['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']], [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_0/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_1/applywarp_v/mapflow/_applywarp_v0/varcope1_warp.nii.gz']]]]
#buffer(copes, varcopes)

def random_combination(iterable):
    import random
    from math import ceil
    out = []
    n = len(iterable)
    r = ceil(n / 2)
    num = range(n)
    
    indices = sorted(random.sample(num, r))
    group = [iterable[i] for i in indices]
    out.append(tuple(group))
    
    out.append(tuple([sub for sub in iterable if sub not in group]))

    return tuple(out)

def make_groups(subjects, split_half, ref):
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
            return group_container
                    
subjects = ['002S6009', '002S6030', '002S6053', '002S6103', '003S4644', '003S6014', '003S6067', '003S6256', '006S0731', '009S6212', '011S6418', '020S6227', '021S4335', '021S6914', '024S6385', '027S6327', '029S4384', '031S0618', '032S6293', '035S6551', '037S4028', '037S6115', '041S6354', '070S4856', '070S6386', '094S6269', '094S6278', '098S6343', '098S6362', '098S6734', '100S6164', '100S6493', '116S0382', '116S6439', '127S4198', '127S4604', '127S6348', '130S6035', '130S6043', '130S6161', '130S6372', '135S4598', '141S6872', '168S6059', '168S6064', '168S6285', '168S6318', '168S6350', '941S6044', '941S6094']
split_half = True
ref = 10

#make_groups(subjects, split_half, ref)

def group_construction(subjects, split_half=True, average=False):
        import itertools, random
        group_container = []
        if split_half:
            prelim = list(itertools.combinations(subjects, round(len(subjects)/2)))
            pre_len = len(prelim)
            for i in range(pre_len):
                if not (pre_len % 2):
                    if i == (pre_len / 2):
                        break
                    else:
                        group_container.append([list(prelim[i]), list(prelim[-(i+1)])])
                else:
                    missed = [sub for sub in subjects if sub not in list(prelim[i])]
                    group_container.append([missed, list(prelim[i])])
        elif average:
            if type(subjects) == list:
                group_container.append([subjects])
            else:
                group_container.append([[subjects]])
        else:
            #PLACEHOLDER
            group_container.append([subjects, ['Single Group']])
            
        return random.sample(group_container, round(len(group_container)*0.1))


def construction(groups, copes, varcopes):
            import re
            subj_proc = len(copes)
            contrasts = len(copes[0])
            merge_contain_c = []
            merge_contain_v = []
            num_copes = []
            
            for group in groups:
                if group[-1][0] == 'Single Group':
                    for i in range(contrasts):
                        merge_contain_c.append([cope[i] for cope in copes])
                        merge_contain_v.append([varcope[i] for varcope in varcopes])
                        num_copes.append(len(group[0]))
                else:
                    for subset in group:
                        cont_c = []
                        cont_v = []
                        for i in range(contrasts):
                            cont_c.append([cope[i] for cope in copes if re.search('_subject_([0-9]+)', cope[i]).group(1) in subset])
                            cont_v.append([varcope[i] for varcope in varcopes if re.search('_subject_([0-9]+)', varcope[i]).group(1) in subset])
                            num_copes.append(len(subset))
                        
                        merge_contain_c.append(cont_c)
                        merge_contain_v.append(cont_v)       
            
            return merge_contain_c, merge_contain_v, num_copes
        
def t_test(covariate, subjects):
            import pandas as pd
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
                    
                    solo = [labels[i] + ' mean', 'T', [labels[i]], [1]] #-> I THINK THIS MIGHT BE AN F CONTRAST
                    contrast = [[labels[i] + '-' + lab, 'T', [labels[i], lab], [1,-1]] if lab != labels[i] else solo for lab in labels]
                    contrasts.append(contrast)
                    
                #NOTE: FOR THE NON-GROUP COVARIATES THEY ARE ADDED AS IS RIGHT NOW -> NO DEMEANING/ORTHOGONALIZATION
                cov = covariates.drop(groupcat, axis=1)
                cat = categories.drop(groupcat)
                
                for c in cat:
                    labels = cov[c].unique()
                    
                    if type(labels[0]) == str:
                        group = cov.groupby(c)
                        num_groups = len(group.count())
                        group_ids = (group.ngroup() + 1).to_list()
                        encoded = group.ngroup().to_list()
                        EVs[c] = encoded
                    else:
                        EVs[c] = cov[c].to_list()
            else:
                single_group = [1] * len(covariates.index)
                label = 'group_mean'
                EVs[label] = single_group
                group_ids = single_group
                contrasts.append([label, 'T', [label], [1]])
                
                    
            return EVs, contrasts, group_ids

#t_test('/Volumes/NewVolume/super_agers/participants.tsv', subjects)
        

def predict(groups, covariate, bold, mask):
            This_may_need_altering = True
            import re
            from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
            from sklearn.multioutput import MultiOutputClassifier
            import pandas as pd
            import numpy as np
            import nibabel as nib
            covariate = pd.read_table(covariate)
            covariates = covariate.set_index('participant_id')
            categories = covariates.columns
            data = covariates.dtypes
            
            reg = covariates.copy()
            
            for label in categories:
                if data[label] == 'O':
                    reg = reg.drop(label, axis=1)
                else:
                    covariates = covariates.drop(label, axis=1)
                    
            regress = False
            classi = False
            
            if reg.shape[1] < len(categories):
                regress = True
                
            if covariates.shape[1] < len(categories):
                classi = True
                    
            
            #covariates = covariate.columns
            num_cov = len(covariates.columns)
            #MAY NEED TO CHANGE SO THAT ADDITIONAL COLUMNS IN PARTICIPANTS ARE USED AS INPUTS AND ONLY DISEASE IS PREDICTED
            pred = []
            #TEST SPEED ON OLD FORMAT AS WELL
            mask_img = nib.load(mask)
            mask_ind = np.nonzero(mask_img.get_fdata())
            
            t_size = [nib.load(f).shape[-1] for f_c in bold for f in f_c]
            
            size = (int(np.ceil(np.ravel(bold).shape[0])), np.shape(mask_ind)[1]*max(t_size))
            dat = np.zeros(size).astype(np.int16)
            subs = []
            ind = 0
            for file_cont in bold:
                for file in file_cont:
                    sub = re.search('_subject_([0-9S]+)', file).group(1)
                    subs.append(sub)
                    dat[ind, :np.shape(mask_ind)[1]*t_size[ind]] = nib.load(file).get_fdata()[mask_ind].reshape(1, -1)
                    ind += 1
                    
            if num_cov:
                if classi:
                    clf_c = RandomForestClassifier(random_state=2021, n_estimators=50)
                if regress:
                    clf_r = RandomForestRegressor(random_state=2021, n_estimators=50)
                #if num_cov > 1:
                #    clf = MultiOutputClassifier(rfc)
                #else:
                #    clf = rfc
                for group in groups:
                    X1 = group[0]
                    X2 = group[1]
                    if classi:
                        y1 = []
                        y2 = []
                    if regress:
                        y1_r = []
                        y2_r = []
                    
                    ind1 = []
                    ind2 = []
                    
                    index = []
                    for subgroup in group:
                        for s in subgroup:
                            index = 'sub-' + s
                            group_ind = [pos for pos, val in enumerate(subs) if val == s]
                            if s in X1:
                                ind1 += group_ind
                                if classi:
                                    y1.append(covariates.loc[index].to_list()*len(group_ind))
                                if regress:
                                    y1_r.append(reg.loc[index].to_list()*len(group_ind))
                            elif s in X2:
                                ind2 += group_ind
                                if classi:
                                    y2.append(covariates.loc[index].to_list()*len(group_ind))
                                if regress:
                                    y2_r.append(reg.loc[index].to_list()*len(group_ind))
                    
                    if classi:
                        Y1 = np.ravel(np.char.array(y1)).reshape(len(ind1), -1)
                        Y2 = np.ravel(np.char.array(y2)).reshape(len(ind1), -1)
                    if regress:
                        Y1_r = np.ravel(np.array(y1_r)).reshape(len(ind1), -1)
                        Y2_r = np.ravel(np.array(y2_r)).reshape(len(ind1), -1)
                        
                    dat1 = dat[tuple(ind1),...]
                    dat2 = dat[tuple(ind2),...]
                    
                    if classi:
                        clf_c = clf_c.fit(dat1, Y1)
                        prediction = clf_c.predict(dat2)
                        pred.append((prediction == Y2).sum() / (Y2.shape[0]*Y2.shape[0]))
                        
                        clf_c = clf_c.fit(dat2, Y2)
                        prediction = clf_c.predict(dat1)
                        pred.append((prediction == Y1).sum() / (Y1.shape[0]*Y1.shape[0]))
                        
                    if regress:
                        clf_r = clf_r.fit(dat1, Y1_r)
                        prediction = clf_r.predict(dat2)
                        pred.append((prediction == Y2_r).sum() / (Y2.shape[0]*Y2.shape[0]))
                        
                        clf_r = clf_r.fit(dat2, Y2_r)
                        prediction = clf_r.predict(dat1)
                        pred.append((prediction == Y1_r).sum() / (Y1.shape[0]*Y1.shape[0]))
    
    
                    #clf = clf.fit(dat2, Y2)
                    #pred.append(clf.score(dat1, Y1))
            else:
                print("NO PARTICIPANTS FILE")
                
            return pred, np.mean(pred)
        
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
        
def parse_xml(xml, goal, mask, egg, yeg):
    import xml.etree.ElementTree as ET
    import re, os, glob
    search = ''
    fsl = os.getenv('FSLDIR')
    for word in goal.split(' '):
        search +=  word + '.*'
        
    ind = 2 - int(re.search('_([0-9])mm', mask).group(1))
    if xml.lower() in 'harvard-oxford':
        for atlas in glob.glob(fsl + '/data/atlases/HarvardOxford*.xml'):
            tree = ET.parse(atlas)
            root = tree.getroot()
        
            for label in root.iter('label'):
                name = label.text
                if re.search(search, name, re.IGNORECASE):
                    file = fsl + '/data/atlases' + root.findall('.//header/images/imagefile')[ind].text + '.nii.gz'
                    index = label.attrib['index']
                    out_name = name
            
        return file, index, out_name
    
    else:
        print('ERROR')

##def resting_state():
 #   if 
from nipype import Workflow, Node, IdentityInterface, Function, MapNode
from nipype.interfaces.fsl import ApplyWarp, ImageMaths
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces import DataSink
import os

random_combination(['01','02','03','05','06'])

HP = 128
TR = 2.5
outliers = '/Volumes/NewVolume/ds000114-new/working_dir/fmri/preprocess/imageproc/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/art/art.sub-02_ses-test_task-fingerfootlips_bold_roi_mcf_st_flirt_outliers.txt'
smoothed = '/Volumes/NewVolume/ds000114-new/working_dir/fmri/preprocess/imageproc/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/smoothnode/smooth/smooth_su/sub-02_ses-test_task-fingerfootlips_bold_roi_mcf_st_flirt_smooth.nii.gz'
mc_par = '/Volumes/NewVolume/ds000114-new/working_dir/fmri/preprocess/imageproc/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/mc/sub-02_ses-test_task-fingerfootlips_bold_roi_mcf.nii.par'

# =============================================================================
# 
# model = Node(SpecifyModel(input_units='secs', parameter_source='FSL'), name='model')
# model.inputs.outlier_files = outliers
# model.inputs.realignment_parameters = mc_par
# model.inputs.time_repetition = TR
# model.inputs.high_pass_filter_cutoff = HP
# model.inputs.functional_runs = smoothed
# model.inputs.event_files = ['/private/var/folders/mx/mztbckq95hzc7px9341hsc480000gn/T/tmp03abik6m/mean_ts/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf_st_flirt_smooth_warp_ts.txt']
# model.run()
# =============================================================================
from nipype.utils.filemanip import save_json
filename = 'dictionary.json'
data = {'test': [0,1,2,3,4,5],
                   'opinion': 'nic cage is good, actually'}
#os.makedirs(os.getcwd() + '/pipeline_0/test')
#save_json(os.getcwd() + '/pipeline_0/test/'+filename, data)
#dat = Node(DataSink(base_directory=os.getcwd()), name='dat')
#dat.inputs.container = 'pipeline'
#setattr(dat.inputs, 'random_number.@num', ['5'])
#setattr(dat.inputs, 'random_number.@boring', ['6'])
#setattr(dat.inputs, 'slice_timings.@egg', ['1','2','3','4','5','6'])
#dat.run()

def count(val):
            #print(i[0])
            count.ind += 1
            return count.ind
count.ind = 0

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
write.inputs.st = ['/Volumes/NewVolume/ds000114-datasink/working_dir/fmri/preprocess/imageproc/_sessions_test/_subject_03/_bbr_True_bbr_type_signed_bet_frac_0.5_bins_640_cost_mutualinfo_discard_4_dof_f_6_frac_mask_0.3_fwhm_4_interp_spline_iso_4_robust_True_susan_True_warplater_True_wm_thresh_0.5/slicetimer/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf_st.nii.gz', '/Volumes/NewVolume/ds000114-datasink/working_dir/fmri/preprocess/imageproc/_sessions_test/_subject_03/_bbr_True_bbr_type_signed_bet_frac_0.5_bins_640_cost_mutualinfo_discard_4_dof_f_6_frac_mask_0.3_fwhm_8_interp_spline_iso_4_robust_True_susan_True_warplater_True_wm_thresh_0.5/slicetimer/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf_st.nii.gz']
write.inputs.sm = ['/Volumes/NewVolume/ds000114-datasink/working_dir/fmri/preprocess/imageproc/_sessions_test/_subject_03/_bbr_True_bbr_type_signed_bet_frac_0.5_bins_640_cost_mutualinfo_discard_4_dof_f_6_frac_mask_0.3_fwhm_4_interp_spline_iso_4_robust_True_susan_True_warplater_True_wm_thresh_0.5/smoothnode/smooth/smooth_su/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf_st_flirt_smooth.nii.gz', '/Volumes/NewVolume/ds000114-datasink/working_dir/fmri/preprocess/imageproc/_sessions_test/_subject_03/_bbr_True_bbr_type_signed_bet_frac_0.5_bins_640_cost_mutualinfo_discard_4_dof_f_6_frac_mask_0.3_fwhm_8_interp_spline_iso_4_robust_True_susan_True_warplater_True_wm_thresh_0.5/smoothnode/smooth/smooth_su/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf_st_flirt_smooth.nii.gz']
write.inputs.mc = ['/Volumes/NewVolume/ds000114-datasink/working_dir/fmri/preprocess/imageproc/_sessions_test/_subject_03/_bbr_True_bbr_type_signed_bet_frac_0.5_bins_640_cost_mutualinfo_discard_4_dof_f_6_frac_mask_0.3_fwhm_4_interp_spline_iso_4_robust_True_susan_True_warplater_True_wm_thresh_0.5/mc/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf.nii.gz.par', '/Volumes/NewVolume/ds000114-datasink/working_dir/fmri/preprocess/imageproc/_sessions_test/_subject_03/_bbr_True_bbr_type_signed_bet_frac_0.5_bins_640_cost_mutualinfo_discard_4_dof_f_6_frac_mask_0.3_fwhm_8_interp_spline_iso_4_robust_True_susan_True_warplater_True_wm_thresh_0.5/mc/sub-03_ses-test_task-fingerfootlips_bold_roi_mcf.nii.gz.par']
write.inputs.base_dir = os.getcwd()
write.inputs.pipeline_st = 0
write.inputs.task = 'fingerfootlips'

#write.run()

class EGG:
    def __init__(self):
        self.pipeline = 0
    
    def run(self):
        pipeline = self.pipeline
        #def count(val, i=[0]):
        #    #print(i[0])
        #    i[0] += 1
        #    return i[0]
        #THIS IS HOW WE WILL UNPACK DICTIONARIES, AND ITERATIVELY CONNECT THEM
        def choose_func(dic, val):
            if dic:
                printer = """def printer(val):\n\tprint(val)\n\treturn """
                for key in val:
                    printer += 'val["{key}"]'.format(key=key) + ', '
                    
            return printer.rstrip(', ')
        
        def sumnow(a,b,c):
            print(a + b + c)
            return a + b + c
                    
        printer = choose_func(True, {'egg': 'yum', 'in': 'tum', 'too': 'indexx'})
        #a, b, c = printer({'egg': 'yum', 'in': 'tum', 'too': 'indexx'})
        p = Node(Function(input_names='val', output_names=['a', 'b', 'c']), name='p')
        p.inputs.function_str = printer
        #printer = choose_func(False)
        p2 = Node(Function(input_names=['a', 'b', 'c'], function=sumnow), name='p2')
            
        
        A = Workflow(name='EGG')
        #A.base_dir = os.getcwd()
        inputnode = Node(IdentityInterface(fields=['iter']), name='inputnode')
        inputnode.inputs.iter = {'egg': 'yum', 'in': 'tum', 'too': 'indexx'}
        outnode = Node(IdentityInterface(fields='iter'), name='outnode')
        A.connect([#(inputnode, outnode, [(('iter', printer), 'iter')]),
                   (inputnode, p, [('iter', 'val')]),
                   (p, p2, [('a', 'a'),
                            ('b', 'b'),
                            ('c','c')]),
                   ])
        #A.connect(inputnode, ('iter', count), outnode, 'iter')
        A.run()
        for i in range(10):
            count(i)
        
        A=3

egg = EGG()
#egg.run()

tester = {'inputnode': {'fields': ['f1', 'f2', 'f3']}, 'printer': {'input_names': ['a', 'b', 'c']}}

A = {'test': [{'egg': 1, 'carrots': 3}, {'tomat': 0, 'yegg': 2}], 'last': [{'bed': 1, 'led': 3}, {'yak': 5, 'snak': 2}]}

def define_print(keys, ind):
    func = "def printer_"
    func += str(ind)
    func += "("
    reuse = ''
    for key in keys:
        reuse += key + ', '
    reuse = reuse.rstrip(', ')
    func += reuse
    func += '):\n\tprint('
    func += reuse + ')'
    
    return func

check = Workflow('check')
for i, key in enumerate(A):
    keys = [k for dic in A[key] for k in list(dic.keys())]
    vars()['input_' + key] = Node(IdentityInterface(fields=keys), name='input_' + key)
    vars()['print_' + key] = Node(Function(input_names=keys), name='print_' + key)
    vars()['print_' + key].inputs.function_str = define_print(keys, i)
    for dic in A[key]:
        for a in dic:
            setattr(vars()['input_' + key].inputs, a, dic[a])
            check.connect(vars()['input_' + key], a, vars()['print_' + key], a)
            
#check.run()


 
from os.path import join as opj
import os
from nipype.interfaces.fsl.maths import MathsCommand
mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
create_seed = Node(MathsCommand(in_file=mask), name='create_seed')
create_seed.inputs.args = '-mul 0 -add 1 -roi {x} 1 {y} 1 {z} 1 0 1'.format(x=45, y=74, z=51)
seed = create_seed.run().outputs.out_file

make_sphere = Node(MathsCommand(in_file=seed), name='make_sphere')
make_sphere.inputs.args = '-kernel sphere {radius} -fmean'.format(radius=5)
#file = make_sphere.run().outputs.out_file

bold = [['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6009/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6030/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6053/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_002S6103/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S4644/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6014/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6067/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_003S6256/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_006S0731/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/level1/appwarps/_subject_009S6212/_ind_0/applywarp_bold/filtered_func_data_warp.nii.gz']]#[['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_04/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_04/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_05/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_05/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_06/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_06/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_07/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_07/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_08/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_08/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_09/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_09/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114-allsubs/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz']]
#[['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_01/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_02/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_03/applywarp_bold/filtered_func_data_warp.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_retest/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/warpflow/_sessions_test/_task_fingerfootlips/_pipeline_0/_subject_10/applywarp_bold/filtered_func_data_warp.nii.gz']]
covariate = '/Volumes/NewVolume/super_agers/participants.tsv'#'/Users/grahamseasons/fMRI/test_data/ds000114/participants.tsv'
groups = [[['002S6009', '003S4644', '003S6067', '006S0731', '009S6212'], ['002S6030', '002S6053', '002S6103', '003S6014', '003S6256']], [['002S6009', '002S6103', '003S4644', '003S6256', '006S0731'], ['002S6030', '002S6053', '003S6014', '003S6067', '009S6212']], [['002S6030', '002S6053', '003S6014', '006S0731', '009S6212'], ['002S6009', '002S6103', '003S4644', '003S6067', '003S6256']], [['002S6009', '003S4644', '003S6256', '006S0731', '009S6212'], ['002S6030', '002S6053', '002S6103', '003S6014', '003S6067']], [['002S6030', '002S6053', '003S4644', '003S6067', '003S6256'], ['002S6009', '002S6103', '003S6014', '006S0731', '009S6212']], [['002S6009', '002S6053', '003S6014', '003S6256', '009S6212'], ['002S6030', '002S6103', '003S4644', '003S6067', '006S0731']], [['002S6009', '002S6030', '002S6053', '003S4644', '009S6212'], ['002S6103', '003S6014', '003S6067', '003S6256', '006S0731']], [['002S6103', '003S4644', '003S6067', '003S6256', '009S6212'], ['002S6009', '002S6030', '002S6053', '003S6014', '006S0731']], [['002S6053', '003S4644', '003S6256', '006S0731', '009S6212'], ['002S6009', '002S6030', '002S6103', '003S6014', '003S6067']], [['002S6053', '002S6103', '003S4644', '003S6256', '009S6212'], ['002S6009', '002S6030', '003S6014', '003S6067', '006S0731']]]#[[['01', '02', '03', '04', '05'], ['06', '07', '08', '09', '10']], [['01', '02', '03', '04', '06'], ['05', '07', '08', '09', '10']], [['01', '02', '03', '04', '07'], ['05', '06', '08', '09', '10']], [['01', '02', '03', '04', '08'], ['05', '06', '07', '09', '10']], [['01', '02', '03', '04', '09'], ['05', '06', '07', '08', '10']], [['01', '02', '03', '04', '10'], ['05', '06', '07', '08', '09']], [['01', '02', '03', '05', '06'], ['04', '07', '08', '09', '10']], [['01', '02', '03', '05', '07'], ['04', '06', '08', '09', '10']], [['01', '02', '03', '05', '08'], ['04', '06', '07', '09', '10']], [['01', '02', '03', '05', '09'], ['04', '06', '07', '08', '10']], [['01', '02', '03', '05', '10'], ['04', '06', '07', '08', '09']], [['01', '02', '03', '06', '07'], ['04', '05', '08', '09', '10']], [['01', '02', '03', '06', '08'], ['04', '05', '07', '09', '10']], [['01', '02', '03', '06', '09'], ['04', '05', '07', '08', '10']], [['01', '02', '03', '06', '10'], ['04', '05', '07', '08', '09']], [['01', '02', '03', '07', '08'], ['04', '05', '06', '09', '10']], [['01', '02', '03', '07', '09'], ['04', '05', '06', '08', '10']], [['01', '02', '03', '07', '10'], ['04', '05', '06', '08', '09']], [['01', '02', '03', '08', '09'], ['04', '05', '06', '07', '10']], [['01', '02', '03', '08', '10'], ['04', '05', '06', '07', '09']], [['01', '02', '03', '09', '10'], ['04', '05', '06', '07', '08']], [['01', '02', '04', '05', '06'], ['03', '07', '08', '09', '10']], [['01', '02', '04', '05', '07'], ['03', '06', '08', '09', '10']], [['01', '02', '04', '05', '08'], ['03', '06', '07', '09', '10']], [['01', '02', '04', '05', '09'], ['03', '06', '07', '08', '10']], [['01', '02', '04', '05', '10'], ['03', '06', '07', '08', '09']], [['01', '02', '04', '06', '07'], ['03', '05', '08', '09', '10']], [['01', '02', '04', '06', '08'], ['03', '05', '07', '09', '10']], [['01', '02', '04', '06', '09'], ['03', '05', '07', '08', '10']], [['01', '02', '04', '06', '10'], ['03', '05', '07', '08', '09']], [['01', '02', '04', '07', '08'], ['03', '05', '06', '09', '10']], [['01', '02', '04', '07', '09'], ['03', '05', '06', '08', '10']], [['01', '02', '04', '07', '10'], ['03', '05', '06', '08', '09']], [['01', '02', '04', '08', '09'], ['03', '05', '06', '07', '10']], [['01', '02', '04', '08', '10'], ['03', '05', '06', '07', '09']], [['01', '02', '04', '09', '10'], ['03', '05', '06', '07', '08']], [['01', '02', '05', '06', '07'], ['03', '04', '08', '09', '10']], [['01', '02', '05', '06', '08'], ['03', '04', '07', '09', '10']], [['01', '02', '05', '06', '09'], ['03', '04', '07', '08', '10']], [['01', '02', '05', '06', '10'], ['03', '04', '07', '08', '09']], [['01', '02', '05', '07', '08'], ['03', '04', '06', '09', '10']], [['01', '02', '05', '07', '09'], ['03', '04', '06', '08', '10']], [['01', '02', '05', '07', '10'], ['03', '04', '06', '08', '09']], [['01', '02', '05', '08', '09'], ['03', '04', '06', '07', '10']], [['01', '02', '05', '08', '10'], ['03', '04', '06', '07', '09']], [['01', '02', '05', '09', '10'], ['03', '04', '06', '07', '08']], [['01', '02', '06', '07', '08'], ['03', '04', '05', '09', '10']], [['01', '02', '06', '07', '09'], ['03', '04', '05', '08', '10']], [['01', '02', '06', '07', '10'], ['03', '04', '05', '08', '09']], [['01', '02', '06', '08', '09'], ['03', '04', '05', '07', '10']], [['01', '02', '06', '08', '10'], ['03', '04', '05', '07', '09']], [['01', '02', '06', '09', '10'], ['03', '04', '05', '07', '08']], [['01', '02', '07', '08', '09'], ['03', '04', '05', '06', '10']], [['01', '02', '07', '08', '10'], ['03', '04', '05', '06', '09']], [['01', '02', '07', '09', '10'], ['03', '04', '05', '06', '08']], [['01', '02', '08', '09', '10'], ['03', '04', '05', '06', '07']], [['01', '03', '04', '05', '06'], ['02', '07', '08', '09', '10']], [['01', '03', '04', '05', '07'], ['02', '06', '08', '09', '10']], [['01', '03', '04', '05', '08'], ['02', '06', '07', '09', '10']], [['01', '03', '04', '05', '09'], ['02', '06', '07', '08', '10']], [['01', '03', '04', '05', '10'], ['02', '06', '07', '08', '09']], [['01', '03', '04', '06', '07'], ['02', '05', '08', '09', '10']], [['01', '03', '04', '06', '08'], ['02', '05', '07', '09', '10']], [['01', '03', '04', '06', '09'], ['02', '05', '07', '08', '10']], [['01', '03', '04', '06', '10'], ['02', '05', '07', '08', '09']], [['01', '03', '04', '07', '08'], ['02', '05', '06', '09', '10']], [['01', '03', '04', '07', '09'], ['02', '05', '06', '08', '10']], [['01', '03', '04', '07', '10'], ['02', '05', '06', '08', '09']], [['01', '03', '04', '08', '09'], ['02', '05', '06', '07', '10']], [['01', '03', '04', '08', '10'], ['02', '05', '06', '07', '09']], [['01', '03', '04', '09', '10'], ['02', '05', '06', '07', '08']], [['01', '03', '05', '06', '07'], ['02', '04', '08', '09', '10']], [['01', '03', '05', '06', '08'], ['02', '04', '07', '09', '10']], [['01', '03', '05', '06', '09'], ['02', '04', '07', '08', '10']], [['01', '03', '05', '06', '10'], ['02', '04', '07', '08', '09']], [['01', '03', '05', '07', '08'], ['02', '04', '06', '09', '10']], [['01', '03', '05', '07', '09'], ['02', '04', '06', '08', '10']], [['01', '03', '05', '07', '10'], ['02', '04', '06', '08', '09']], [['01', '03', '05', '08', '09'], ['02', '04', '06', '07', '10']], [['01', '03', '05', '08', '10'], ['02', '04', '06', '07', '09']], [['01', '03', '05', '09', '10'], ['02', '04', '06', '07', '08']], [['01', '03', '06', '07', '08'], ['02', '04', '05', '09', '10']], [['01', '03', '06', '07', '09'], ['02', '04', '05', '08', '10']], [['01', '03', '06', '07', '10'], ['02', '04', '05', '08', '09']], [['01', '03', '06', '08', '09'], ['02', '04', '05', '07', '10']], [['01', '03', '06', '08', '10'], ['02', '04', '05', '07', '09']], [['01', '03', '06', '09', '10'], ['02', '04', '05', '07', '08']], [['01', '03', '07', '08', '09'], ['02', '04', '05', '06', '10']], [['01', '03', '07', '08', '10'], ['02', '04', '05', '06', '09']], [['01', '03', '07', '09', '10'], ['02', '04', '05', '06', '08']], [['01', '03', '08', '09', '10'], ['02', '04', '05', '06', '07']], [['01', '04', '05', '06', '07'], ['02', '03', '08', '09', '10']], [['01', '04', '05', '06', '08'], ['02', '03', '07', '09', '10']], [['01', '04', '05', '06', '09'], ['02', '03', '07', '08', '10']], [['01', '04', '05', '06', '10'], ['02', '03', '07', '08', '09']], [['01', '04', '05', '07', '08'], ['02', '03', '06', '09', '10']], [['01', '04', '05', '07', '09'], ['02', '03', '06', '08', '10']], [['01', '04', '05', '07', '10'], ['02', '03', '06', '08', '09']], [['01', '04', '05', '08', '09'], ['02', '03', '06', '07', '10']], [['01', '04', '05', '08', '10'], ['02', '03', '06', '07', '09']], [['01', '04', '05', '09', '10'], ['02', '03', '06', '07', '08']], [['01', '04', '06', '07', '08'], ['02', '03', '05', '09', '10']], [['01', '04', '06', '07', '09'], ['02', '03', '05', '08', '10']], [['01', '04', '06', '07', '10'], ['02', '03', '05', '08', '09']], [['01', '04', '06', '08', '09'], ['02', '03', '05', '07', '10']], [['01', '04', '06', '08', '10'], ['02', '03', '05', '07', '09']], [['01', '04', '06', '09', '10'], ['02', '03', '05', '07', '08']], [['01', '04', '07', '08', '09'], ['02', '03', '05', '06', '10']], [['01', '04', '07', '08', '10'], ['02', '03', '05', '06', '09']], [['01', '04', '07', '09', '10'], ['02', '03', '05', '06', '08']], [['01', '04', '08', '09', '10'], ['02', '03', '05', '06', '07']], [['01', '05', '06', '07', '08'], ['02', '03', '04', '09', '10']], [['01', '05', '06', '07', '09'], ['02', '03', '04', '08', '10']], [['01', '05', '06', '07', '10'], ['02', '03', '04', '08', '09']], [['01', '05', '06', '08', '09'], ['02', '03', '04', '07', '10']], [['01', '05', '06', '08', '10'], ['02', '03', '04', '07', '09']], [['01', '05', '06', '09', '10'], ['02', '03', '04', '07', '08']], [['01', '05', '07', '08', '09'], ['02', '03', '04', '06', '10']], [['01', '05', '07', '08', '10'], ['02', '03', '04', '06', '09']], [['01', '05', '07', '09', '10'], ['02', '03', '04', '06', '08']], [['01', '05', '08', '09', '10'], ['02', '03', '04', '06', '07']], [['01', '06', '07', '08', '09'], ['02', '03', '04', '05', '10']], [['01', '06', '07', '08', '10'], ['02', '03', '04', '05', '09']], [['01', '06', '07', '09', '10'], ['02', '03', '04', '05', '08']], [['01', '06', '08', '09', '10'], ['02', '03', '04', '05', '07']], [['01', '07', '08', '09', '10'], ['02', '03', '04', '05', '06']]]#[[['01', '02'], ['03', '10']], [['01', '03'], ['02', '10']], [['01', '10'], ['02', '03']]]
mask = opj(os.getenv('FSLDIR'), 'data/linearMNI/MNI152lin_T1_2mm_brain.nii.gz')
xml = 'harvard'#opj(os.getenv('FSLDIR'), 'data/atlases/HarvardOxford-Cortical.xml')
predict(groups, covariate, bold, mask)
goal = 'cing post'#'inf temp temp'
dicti = {'egg': 'yum', 'yeg': 'yuck'}
#parse_xml(xml, goal, mask, **dicti)

stats = [['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo0/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo1/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo2/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo3/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo4/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo5/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo6/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo7/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo8/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo9/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo10/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo11/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo12/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo13/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo14/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo15/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo16/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo17/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo18/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/sup_test/working_dir/fmri/l2analysis/l2/_ind_0/flameo/mapflow/_flameo19/stats/zstat1.nii.gz']]#[['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo0/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo1/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo2/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo3/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo4/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo5/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo6/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo7/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo8/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo9/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo10/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo11/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo12/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo13/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo14/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo15/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo16/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo17/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo18/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo19/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo20/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo21/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo22/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo23/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo24/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo25/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo26/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo27/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo28/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo29/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo30/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo31/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo32/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo33/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo34/stats/zstat1.nii.gz'], ['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo35/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo36/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo37/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo38/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo39/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo40/stats/zstat1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis/l2/_task_fingerfootlips/_pipeline_0/flameo/mapflow/_flameo41/stats/zstat1.nii.gz']]
#compare(stats, mask)    
#copes = [[['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo0/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo1/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo2/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo3/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo4/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo5/stats/cope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo6/stats/cope1.nii.gz']]]
#groups = [[['01', '02'], ['03', '10']], [['01', '03'], ['02', '10']], [['01', '10'], ['02', '03']]]
#varcopes = [[['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_01/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_02/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_03/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']], [['/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo0/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo1/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo2/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo3/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo4/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo5/stats/varcope1.nii.gz', '/Volumes/NewVolume/ds000114/working_dir/fmri/l2analysis_fe/l2/_task_fingerfootlips/_pipeline_0/_subject_10/_subgroup_0/flameo/mapflow/_flameo6/stats/varcope1.nii.gz']]]
#t_test('/Users/grahamseasons/fMRI/test_data/ds000114/participants.tsv', ['01', '02', '03', '04'])
#construction(groups, copes, varcopes)
#group_construction(['01','02','03','04'])