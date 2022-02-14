#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 12:50:07 2021

@author: grahamseasons
"""
def check4brains(data_dir, T1w):
    from nipype import SelectFiles
    from niworkflows.anat.ants import init_brain_extraction_wf
    from nipype import DataSink, Node
    import re, os, glob
    base_dir = os.getcwd()
    
    sID = re.search('sub-([0-9S]+)', T1w)
    if sID:
        sID = sID.group(1)
    else:
        raise NameError('Subject ID not found. ID should be composed of digits, with one letter (S) allowed')
    
    templates = {'bet': data_dir + '/brain_extracted/_subject_{sID}/*masked.nii.gz'.format(sID=sID)}
    try:
        selectfiles = SelectFiles(templates, sort_filelist=True, force_lists=['bet']).run().outputs.bet
    except:
        selectfiles = []
        
    if selectfiles:
        return selectfiles
    
    data = Node(DataSink(parameterization=False), name='data', base_dir=base_dir)
    
    wf = init_brain_extraction_wf()
    wf.base_dir = os.getcwd()
    wf.inputs.inputnode.in_files = T1w
    wf.connect(wf.get_node('outputnode'), 'out_file', data, data_dir + '/brain_extracted/_subject_{sID}/'.format(sID=sID))
    wf.run()
    
    return glob.glob(data_dir + '/brain_extracted/_subject_{sID}/*masked.nii*'.format(sID=sID))[0]

def smooth(warped, mask):
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
                
        smooth.run()
        smoothed = glob.glob(os.getcwd() + '/smooth/smooth_su/*.nii*')[0]
    else:
        smooth_iso = Node(IsotropicSmooth(fwhm=fwhm, output_type='NIFTI'), name='smooth_iso')
        smooth.connect(inputnode, 'in_file', smooth_iso, 'in_file')
        smooth.connect(smooth_iso, 'out_file', outnode, 'smoothed')
        
        smooth.run()
        smoothed = glob.glob(os.getcwd() + '/smooth/smooth_iso/*.nii*')[0]
        
    return smoothed, glob.glob(os.getcwd() + '/smooth/**', recursive=True)
                        
def registration(T1w, mask, start_img, corrected_img, bet, wm_file):
    from nipype.interfaces.fsl import FLIRT
    from nipype.interfaces.base import Undefined
    from nipype import IdentityInterface
    from nipype import Workflow, Node
    from os.path import join as opj
    import os, glob, re
    
    reg = Workflow(name='reg')
    reg.base_dir = os.getcwd()
    
    bbr = vars().get('bbr', True)
    concatenate = vars().get('concatenate', True)
    #CHECK TO VERIFY THAT DOF IS CALLED DOG
    regpre = Node(FLIRT(in_file=start_img, dof=6, reference=bet, output_type='NIFTI_GZ'), name='regpre')
    
    applywarp = Node(FLIRT(in_file=corrected_img, reference=bet, output_type='NIFTI_GZ'), name='applywarp')
    
    outnode =  Node(IdentityInterface(fields=['warped','out_mat', 'mask']), name='outnode')
    
    if bbr:
        regbbr = Node(FLIRT(cost='bbr', reference=T1w, in_file=start_img, output_type='NIFTI_GZ',
                        schedule=opj(os.getenv('FSLDIR'), 'etc/flirtsch/bbr.sch')), name='regbbr')
        regbbr.inputs.wm_seg = wm_file
        
        reg.connect([(regpre, regbbr, [('out_matrix_file', 'in_matrix_file')]),
                     (regbbr, outnode, [('out_matrix_file', 'out_mat')]),
                     ])
        
        if not concatenate:
            reg.connect([(regbbr, applywarp, [('out_matrix_file', 'in_matrix_file')])])
            
        node_reg = 'regbbr'
    else:
        reg.connect([(regpre, outnode, [('out_matrix_file', 'out_mat')]),
                     ])
        
        if not concatenate:
            reg.connect([(regpre, applywarp, [('out_matrix_file', 'in_matrix_file')])])
            
        node_reg = 'regpre'

    node_warp = 'applywarp'
        
    reg.run()
    
    out_mat = glob.glob(os.getcwd() + '/reg/' + node_reg + '/*.mat')[0]
    if concatenate:
        warped = corrected_img
    else:
        warped = glob.glob(os.getcwd() + '/reg/' + node_warp + '/*.nii*')[0]
        
    return out_mat, warped, glob.glob(os.getcwd() + '/reg/**', recursive=True)


def regress(unsmoothed, mc_par, segmentations, mask, rest):
    from nipype import Node
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.fsl import ImageMeants, Threshold, FLIRT, GLM, ImageStats, ImageMaths
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
    from nipype.interfaces.fsl import ApplyWarp, ConvertWarp, Threshold, ConvertXFM, FLIRT
    import re
    from nipype.interfaces.base import Undefined
    if not vars().get('warplater', True):
        
        for i, file in enumerate(segmentations):
            segmentations[i] = ApplyWarp(in_file=file, ref_file=mniMask, field_file=warp, interp='nn').run().outputs.out_file
        
        if vars().get('concatenate', True):
            warp = ConvertWarp(reference=mniMask, premat=out_mat, warp1=warp).run().outputs.out_file
            
        warped = ApplyWarp(in_file=warped, ref_file=mniMask, field_file=warp)
        
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
    
    
    
    
    
    