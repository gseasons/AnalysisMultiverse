#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 12:50:07 2021

@author: grahamseasons
"""
def betwrite(data_dir, bet_):
    """Writes brain extracted image to permanent location"""
    from nipype import DataSink, Node
    import re, os
    base_dir = os.getcwd()
    data = Node(DataSink(parameterization=False), name='data', base_dir=base_dir)
    sID = re.search('sub-([0-9S]+)', bet_)
    setattr(data.inputs, data_dir + '/brain_extracted/_subject_{sID}/'.format(sID=sID), [bet_])
    data.run()

def check4brains(data_dir, T1w):
    """Checks to see if there are saved brain extractions"""
    from nipype import SelectFiles
    import re
    
    sID = re.search('sub-([0-9S]+)', T1w)
    if sID:
        sID = sID.group(1)
    else:
        raise NameError('Subject ID not found. ID should be composed of digits, with one letter (S) allowed')
    
    templates = {'bet': data_dir + '/brain_extracted/_subject_{sID}/*masked.nii.gz'.format(sID=sID)}
    try:
        selectfiles = SelectFiles(templates, sort_filelist=True, force_lists=['bet']).run().outputs.bet
    except:
        raise FileNotFoundError('Missing brain extracted image for subject {0}'.format(sID))
        
    return selectfiles

def smooth(warped, mask):
    """Implements both FSL SUSAN smoothing, and Isometric smoothing"""
    #WORKFLOW ADAPTED FROM: https://nipype.readthedocs.io/en/latest/users/examples/fmri_fsl.html
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.fsl import (ImageMaths, ImageStats, SUSAN, IsotropicSmooth)
    import os, glob, re
    from preprocessing.functions import get_bright_thresh, getthreshop
    
    susan = vars().get('susan', True)
    fwhm = vars().get('fwhm', 6)
    
    if susan:
        maskfunc = ImageMaths(suffix='_bet', op_string='-mas', in_file2=mask, in_file=warped).run().outputs.out_file
        getthresh = ImageStats(op_string='-p 2 -p 98', in_file=maskfunc).run().outputs.out_stat
        getthreshop_ = getthreshop(getthresh)
        threshold = ImageMaths(out_data_type='char', suffix='_thresh', in_file=maskfunc, op_string=getthreshop_).run().outputs.out_file
        
        medianval = ImageStats(op_string='-k %s -p 50', in_file=warped, mask_file=threshold).run().outputs.out_stat
        
        brightthresh = get_bright_thresh(medianval)
        
        smoothed = SUSAN(fwhm=fwhm, brightness_threshold=brightthresh, in_file=warped).run().outputs.smoothed_file
    else:
        smoothed = IsotropicSmooth(fwhm=fwhm, in_file=warped, output_type='NIFTI').run().outputs.out_file
        
    return smoothed
                        
def registration(T1w, start_img, corrected_img, bet, wm_file):
    """Coregsitration of functional image to T1w image"""
    #NOTE: MASK NOT USED -> COULD REMOVE
    from nipype.interfaces.fsl import FLIRT
    from nipype.interfaces.base import Undefined
    from nipype import Node
    from os.path import join as opj
    import os, re
    base_dir = os.getcwd()
    
    bbr = vars().get('bbr', True)
    concatenate = vars().get('concatenate', True)
    
    regpre = Node(FLIRT(in_file=start_img, dof=6, reference=bet, output_type='NIFTI_GZ'), name='regpre', base_dir=base_dir)
    applywarp = Node(FLIRT(in_file=corrected_img, reference=bet, output_type='NIFTI_GZ'), name='applywarp', base_dir=base_dir)
    
    if bbr:
        regbbr = Node(FLIRT(cost='bbr', reference=T1w, in_file=start_img, output_type='NIFTI_GZ',
                        schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch')), name='regbbr', base_dir=base_dir)
        regbbr.inputs.wm_seg = wm_file
        
        regpre = regpre.run().outputs.out_matrix_file
        
        regbbr.inputs.in_matrix_file = regpre
        outmat = regbbr.run().outputs.out_matrix_file
    else:
        
        outmat = regpre.run().outputs.out_matrix_file
    
    if not concatenate:
        applywarp.inputs.in_matrix_file = outmat
        warped = applywarp.run().outputs.out_file
    else:
        warped = corrected_img
        
    return outmat, warped


def regress(unsmoothed, mc_par, segmentations, mask, rest):
    """Signal regression (WM, CSF, Global, motion)"""
    from nipype import Node
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.fsl import ImageMeants, Threshold, FLIRT, GLM, ImageStats, ImageMaths
    from nipype.interfaces.fsl.maths import MeanImage
    from nipype.interfaces.base import CommandLine
    import numpy as np
    import nibabel as nib
    import os, re
    base_dir = os.getcwd()
    #For data driven analysis, regression is constant according to paper it is based on
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
            
    meanglob = Node(ImageMeants(in_file=unsmoothed, mask=mask), name='meanglob', base_dir=base_dir)
    glob = np.loadtxt(meanglob.run().outputs.out_file).reshape(-1, 1)
    reho = np.hstack((reho, glob))
        
    if GLOBAL:
        params = np.hstack((params, glob))
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
    
        regressed = glm.run().outputs.out_res
        add = abs(ImageStats(in_file=regressed, op_string='-R').run().outputs.out_stat[0])
        unsmoothed = ImageMaths(in_file=regressed, args='-add '+str(add)).run().outputs.out_file
        
    if np.any(reho) and not np.array_equal(reho, params):
        out_name = re.search('/([A-Za-z0-9_-]+).nii', forreho).group(1) + '_reho_regressed.nii.gz'
        glm = GLM(design=name_reho + '.mat', in_file=forreho, out_res_name=out_name)
    
        regressed = glm.run().outputs.out_res
        add = abs(ImageStats(in_file=regressed, op_string='-R').run().outputs.out_stat[0])
        forreho = ImageMaths(in_file=regressed, args='-add '+str(add)).run().outputs.out_file
        
    return unsmoothed, forreho


def mni(mniMask, brain, brainmask, warp, warped, segmentations, out_mat, start_img):
    """Converts to MNI space or the inverse if conversion takes place after l1 analysis"""
    from nipype.interfaces.fsl import ApplyWarp, ConvertWarp, Threshold, ConvertXFM, FLIRT, Merge
    from nipype.interfaces.fsl.maths import MaxImage
    import re
    from nipype.interfaces.base import Undefined
    
    if not vars().get('warplater', True):
        for i, file in enumerate(segmentations):
            segmentations[i] = ApplyWarp(in_file=file, ref_file=mniMask, field_file=warp, interp='nn').run().outputs.out_file
        
        if vars().get('concatenate', True):
            warp = ConvertWarp(reference=mniMask, premat=out_mat, warp1=warp).run().outputs.out_file
            
        warped = ApplyWarp(in_file=warped, ref_file=mniMask, field_file=warp)
        start_warped = ApplyWarp(in_file=start_img, ref_file=mniMask, field_file=warp)
        
        warped = warped.run().outputs.out_file
        startimg = start_warped.run().outputs.out_file
        if brainmask == '':
            brainmask = Merge(in_files=segmentations, dimension='t').run().outputs.merged_file
            brainmask = MaxImage(in_file=brainmask, dimension='T').run().outputs.out_file
        else:
            brainmask = ApplyWarp(in_file=brainmask, ref_file=mniMask, field_file=warp, interp='nn').run().outputs.out_file
    else:
        startimg = start_img
        if not vars().get('concatenate', False):
            brainmask = Threshold(in_file=brain, thresh=0, args='-bin').run().outputs.out_file
        else:
            warp = ConvertWarp(reference=mniMask, premat=out_mat, warp1=warp).run().outputs.out_file
            segform = ConvertXFM(in_file=out_mat, invert_xfm=True).run().outputs.out_file
            for i, file in enumerate(segmentations):
                segmentations[i] = FLIRT(in_file=file, reference=warped, interp='nearestneighbour', apply_xfm=True, in_matrix_file=segform).run().outputs.out_file
            
            if brainmask == '':
                brainmask = Merge(in_files=segmentations, dimension='t').run().outputs.merged_file
                brainmask = MaxImage(in_file=brainmask, dimension='T').run().outputs.out_file
            
    return warped, brainmask, segmentations, warp, startimg
    
    
    
    
    
    