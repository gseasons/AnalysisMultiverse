#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:06:49 2021

@author: grahamseasons
"""
from nipype import Node, Workflow, MapNode, IdentityInterface, Function, SelectFiles
from nipype.interfaces.fsl import FLIRT, FNIRT, ApplyWarp, DilateImage, UnaryMaths
import os

class spatial_normalization:
    def __init__(self):
        self.dec = True
        
    def construct(self):
        if self.dec:
            self.dec = False
            return self.get_warps()
        else:
            self.dec = True
            return self.apply_warps()
    
    def get_warps(self, flow=''):
        from updated.normalization.functions import invert
        inputnode = Node(IdentityInterface(fields=['brain', 'boldmask', 'concatenate', 'ref_file']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['warp', 'invwarp']), name='outnode')
        
        prelim = Node(FLIRT(dof=12, output_type='NIFTI_GZ'), name='prelim')
        warp = Node(FNIRT(field_file=True, config_file='T1_2_MNI152_2mm'), name='warp')
        invwarp = Node(Function(input_names=['warp', 'brain', 'brainmask', 'warplater', 'coregmat', 'concatenate'],
                                output_names=['invwarp'], function=invert), name='invwarp')
        
        dilatebrain = Node(DilateImage(operation='max'), name='dilatebrain')
        dilateref = Node(DilateImage(operation='max'), name='dilateref')
        binarizebrain = Node(UnaryMaths(operation='bin'), name='binarizebrain')
        binarizeref = Node(UnaryMaths(operation='bin'), name='binarizeref')
        
        if flow:
            genwarps = flow
        else:
            genwarps = Workflow('genwarps')
            genwarps.base_dir = os.getcwd()
            genwarps.connect([(inputnode, prelim, [('brain', 'in_file')]),
                              (inputnode, prelim, [('ref_file', 'reference')]),
                              (inputnode, warp, [('brain', 'in_file')]),
                              (inputnode, warp, [('ref_file', 'ref_file')]),
                              (inputnode, invwarp, [('warppostfeat', 'warplater')]),
                              (inputnode, invwarp, [('boldmask', 'brainmask')]),
                              (inputnode, invwarp, [('no_resample', 'no_resample')]),
                              (inputnode, invwarp, [('brain', 'brain')]),
                              (warp, invwarp, [('field_file', 'warp')]),
                              (warp, outnode, [('field_file', 'warp')]),
                              (invwarp, outnode, [('invwarp', 'invwarp')]),
                              ])
        
        genwarps.connect([(prelim, warp, [('out_matrix_file', 'affine_file')]),
                          (dilatebrain, binarizebrain, [('out_file', 'in_file')]),
                          (dilateref, binarizeref, [('out_file', 'in_file')]),
                          (binarizebrain, warp, [('out_file', 'inmask_file')]),
                          (binarizeref, warp, [('out_file', 'refmask_file')]),
                          (warp, invwarp, [('field_file', 'warp')]),
                          ])
        return genwarps
        
    def apply_warps(self, flow=''):
        from updated.normalization.functions import identity, ret_files, check_bold
        inputnode = Node(IdentityInterface(fields=['feat_dir', 'warp_file', 'ref_file', 'needwarp']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['feat_dir', 'warp_file', 'ref_file', 'needwarp']), name='outnode')
        
        templates = {'cope': 'stats/cope*.nii.gz',
                     'varcope': 'stats/varcope*.nii.gz',
                     'bold': 'filtered_func_data.nii.gz'}
        
        selectfiles = Node(SelectFiles(templates, sort_filelist=True), name='selectfiles')
        
        ident = MapNode(Function(input_names=['cope', 'varcope', 'bold', 'needwarp'],
                              output_names=['cope', 'varcope', 'bold'], function=identity), name='ident', iterfield=['cope', 'varcope'])
            
        ret = MapNode(Function(input_names=['cope_orig', 'varcope_orig', 'bold_orig', 'cope_warp', 
                                         'varcope_warp', 'bold_warp', 'needwarp'],
                               output_names=['cope', 'varcope', 'bold'], function=ret_files), name='ret', iterfield=['cope_orig', 'varcope_orig', 'cope_warp', 'varcope_warp'])

        applywarpcopes = MapNode(ApplyWarp(), name='applywarpcopes', iterfield=['in_file'])
        applywarpvarcopes = MapNode(ApplyWarp(), name='applywarpvarcopes', iterfield=['in_file'])
        applywarpbold = Node(ApplyWarp(), name='applywarpbold')
        
        if flow:
            appwarps = flow
        else:
            appwarps = Workflow('genwarps')
            appwarps.base_dir = os.getcwd()
            appwarps.connect([(inputnode, applywarpcopes, [('ref_file', 'ref_file')]),
                              (inputnode, applywarpvarcopes, [('ref_file', 'ref_file')]),
                              (inputnode, applywarpbold, [('ref_file', 'ref_file')]),
                              (inputnode, applywarpcopes, [('warp_file', 'field_file')]),
                              (inputnode, applywarpvarcopes, [('warp_file', 'field_file')]),
                              (inputnode, applywarpbold, [('warp_file', 'field_file')]),
                              (inputnode, selectfiles, [('feat_dir', 'feat_dir')]),
                              (inputnode, ret, [('needwarp', 'needwarp')]),
                              (inputnode, ident, [('needwarp', 'needwarp')]),
                              (ret, outnode, [('cope', 'cope')]),
                              (ret, outnode, [('varcope', 'varcope')]),
                              (ret, outnode, [(('bold', check_bold), 'bold')]),
                              ])
            
        appwarps.connect([(selectfiles, ident, [('cope', 'cope')]),
                          (selectfiles, ident, [('varcope', 'varcope')]),
                          (selectfiles, ident, [('bold', 'bold')]),
                          (ident, applywarpcopes, [('cope', 'in_file')]),
                          (ident, applywarpvarcopes, [('varcope', 'in_file')]),
                          (ident, applywarpbold, [(('bold', check_bold), 'in_file')]),
                          (applywarpcopes, ret, [('out_file', 'cope_warp')]),
                          (applywarpvarcopes, ret, [('out_file', 'varcope_warp')]),
                          (applywarpbold, ret, [('out_file', 'bold_warp')]),
                          (selectfiles, ret, [('cope', 'cope_orig')]),
                          (selectfiles, ret, [('varcope', 'varcope_orig')]),
                          (selectfiles, ret, [('bold', 'bold_orig')]),
                          ])
        return appwarps