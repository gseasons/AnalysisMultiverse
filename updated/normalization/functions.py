#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:31:37 2021

@author: grahamseasons
"""
def invert(warp, brain, warplater):
    from nipype import Node
    from nipype import InvWarp
    if warplater:
        invwarp = Node(InvWarp(warp=warp, ref=brain), name='invwarp')
        invwarp = invwarp.run().outputs.inverse_warp
    else:
        invwarp = ''
    return invwarp

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

def ret_files(cope_orig, varcope_orig, bold_orig, cope_warp, varcope_warp, bold_warp, needwarp):
    if not needwarp:
        return cope_orig, varcope_orig, bold_orig
    else:
        return cope_warp, varcope_warp, bold_warp
    
def check_bold(bolds):
    if type(bolds) == list:
        return bolds[0]
    else:
        return bolds