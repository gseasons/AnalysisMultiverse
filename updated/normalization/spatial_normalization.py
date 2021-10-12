#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:06:49 2021

@author: grahamseasons
"""
from nipype import Node, Workflow, MapNode, IdentityInterface, Function, SelectFiles, JoinNode
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
        inputnode = Node(IdentityInterface(fields=['brain', 'ref_file']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['warp']), name='outnode')
        
        prelim = Node(FLIRT(dof=12, output_type='NIFTI_GZ'), name='prelim')
        warp = Node(FNIRT(field_file=True, config_file='T1_2_MNI152_2mm'), name='warp')
        
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
                              (warp, outnode, [('field_file', 'warp')]),
                              ])
        
        genwarps.connect([(prelim, warp, [('out_matrix_file', 'affine_file')]),
                          (dilatebrain, binarizebrain, [('out_file', 'in_file')]),
                          (dilateref, binarizeref, [('out_file', 'in_file')]),
                          (binarizebrain, warp, [('out_file', 'inmask_file')]),
                          (binarizeref, warp, [('out_file', 'refmask_file')]),
                          ])
        if flow:
            return genwarps
        
    def apply_warps(self, apply_var):
        appwarps = Workflow('appwarps')
        appwarps.base_dir = os.getcwd()
        
        apply_var = apply_var['apply']
        
        #inputnode = Node(IdentityInterface(fields=['feat_dir', 'warp_file', 'ref_file', 'needwarp']), name='inputnode')
        
        rng = len(apply_var['PIPELINE'][0][1])
        
        fields=['feat_dir', 'warp_file', 'ref_file', 'ind']
        
        fields += [key[0] for key in apply_var['WARP'] if key != 'PIPELINE']
        
        inputnode = Node(IdentityInterface(fields=fields), name='inputnode_app')
        inputnode.iterables = apply_var['ind']
        inputnode.synchronize = True
        
        for tup in apply_var['WARP']:
            setattr(inputnode.inputs, tup[0], tup[1])
        
        if rng > 1:
            outnode = JoinNode(IdentityInterface(fields=['cope', 'varcope', 'bold']), name='outnode', joinsource='inputnode_app', joinfield=['cope', 'varcope', 'bold'])
        else:    
            outnode = Node(IdentityInterface(fields=['cope', 'varcope', 'bold']), name='outnode')
        
        templates = {'cope': 'stats/cope*.nii.gz',
                     'varcope': 'stats/varcope*.nii.gz',
                     'bold': 'filtered_func_data.nii.gz'}
        
        selectfiles = Node(SelectFiles(templates, sort_filelist=True), name='selectfiles')
        
        def identity(cope, varcope, bold, needwarp):
            from nipype.interfaces.fsl import ImageMaths
            from nipype import Node
            
            if not needwarp:
                clear = Node(ImageMaths(op_string='-mul 0 -bin'), name='clear')
                clear.inputs.in_file = cope
                cope = clear.run().outputs.out_file
                varcope = cope
                bold = cope
            
            return cope, varcope, bold
        
        ident = MapNode(Function(input_names=['cope', 'varcope', 'bold', 'needwarp'],
                              output_names=['cope', 'varcope', 'bold'], function=identity), name='ident', iterfield=['cope', 'varcope'])
        
        def ret_files(cope_orig, varcope_orig, bold_orig, cope_warp, varcope_warp, bold_warp, needwarp):
            if not needwarp:
                return cope_orig, varcope_orig, bold_orig
            else:
                return cope_warp, varcope_warp, bold_warp
            
        ret = MapNode(Function(input_names=['cope_orig', 'varcope_orig', 'bold_orig', 'cope_warp', 
                                         'varcope_warp', 'bold_warp', 'needwarp'],
                            output_names=['cope', 'varcope', 'bold'], function=ret_files), name='ret', iterfield=['cope_orig', 'varcope_orig', 'cope_warp', 'varcope_warp'])
        
        def check_bold(bolds):
            if type(bolds) == list:
                return bolds[0]
            else:
                return bolds
            
        def buffer(feat_dir, needwarp, i, rng):
            if rng > 1:
                return feat_dir[i], needwarp[i]
            else:
                return feat_dir, needwarp
        
        buff = Node(Function(input_names=['feat_dir', 'needwarp', 'i', 'rng'],
                             output_names=['feat_dir', 'needwarp'], function=buffer), name='buff')
        buff.inputs.rng = rng
        
        def single_bold(bold):
            return bold[0]
        
        applywarp_c = MapNode(ApplyWarp(), name='applywarp_c', iterfield=['in_file'])
        applywarp_v = MapNode(ApplyWarp(), name='applywarp_v', iterfield=['in_file'])
        applywarp_bold = Node(ApplyWarp(), name='applywarp_bold')
        
        appwarps.connect([(inputnode, applywarp_c, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_v, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_bold, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_c, [('warp_file', 'field_file')]),
                          (inputnode, applywarp_v, [('warp_file', 'field_file')]),
                          (inputnode, applywarp_bold, [('warp_file', 'field_file')]),
                          (inputnode, buff, [('ind', 'i')]),
                          (inputnode, buff, [('feat_dir', 'feat_dir')]),
                          (buff, selectfiles, [('feat_dir', 'base_directory')]),
                          (inputnode, buff, [('needwarp', 'needwarp')]),
                          (buff, ret, [('needwarp', 'needwarp')]),
                          (selectfiles, ident, [('cope', 'cope')]),
                          (selectfiles, ident, [('varcope', 'varcope')]),
                          (selectfiles, ident, [('bold', 'bold')]),
                          (buff, ident, [('needwarp', 'needwarp')]),
                          (ident, applywarp_c, [('cope', 'in_file')]),
                          (ident, applywarp_v, [('varcope', 'in_file')]),
                          (ident, applywarp_bold, [(('bold', check_bold), 'in_file')]),
                          (applywarp_c, ret, [('out_file', 'cope_warp')]),
                          (applywarp_v, ret, [('out_file', 'varcope_warp')]),
                          (applywarp_bold, ret, [('out_file', 'bold_warp')]),
                          (selectfiles, ret, [('cope', 'cope_orig')]),
                          (selectfiles, ret, [('varcope', 'varcope_orig')]),
                          (selectfiles, ret, [('bold', 'bold_orig')]),
                          (ret, outnode, [('cope', 'cope')]),
                          (ret, outnode, [('varcope', 'varcope')]),
                          (ret, outnode, [(('bold', single_bold), 'bold')]),
                          ])
        
        return appwarps