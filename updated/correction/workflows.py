#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 12:47:09 2021

@author: grahamseasons
"""
def FDR(zstat, mask, cor):
    from updated.correction.functions import fdr
    from nipype import MapNode, Function
    from nipype.interfaces.fsl import ImageMaths, Threshold, BinaryMaths
    import os
    
    cor['q'] = cor.pop('p')

    pos_thresh = MapNode(Threshold(direction='below'), name='pos_thresh', iterfield='in_file', nested=True)
    pos_thresh.inputs.thresh = 0
    pos_thresh.inputs.in_file = zstat
    pthreshed = pos_thresh.run().outputs.out_file
    
    neg_thresh = MapNode(Threshold(direction='above'), name='neg_thresh', iterfield='in_file', nested=True)
    neg_thresh.inputs.thresh = 0
    neg_thresh.inputs.in_file = zstat
    nthreshed = neg_thresh.run().outputs.out_file
    
    p_file_pos = MapNode(ImageMaths(op_string='-ztop', suffix='_pval'), 
                     name='p_file_pos', iterfield=['in_file'], nested=True)
    p_file_pos.inputs.in_file = pthreshed
    p_file_pos = p_file_pos.run().outputs.out_file
    
    p_file_neg = MapNode(ImageMaths(op_string='-ztop -mul -1 -add 1', suffix='_pval'), 
                     name='p_file_neg', iterfield=['in_file'], nested=True)
    p_file_neg.inputs.in_file = nthreshed
    p_file_neg = p_file_neg.run().outputs.out_file
    
    form_pos = MapNode(Function(input_names=['p_im', 'mask', 'q'],
                                output_names=['form_str'], function=fdr), 
                       name='form_pos', iterfield='p_im', nested=True)
    form_pos.inputs.p_im = p_file_pos
    form_pos.inputs.mask = mask
    form_pos.inputs.q = cor['q']
    pos = form_pos.run().outputs.form_str
    
    form_neg = MapNode(Function(input_names=['p_im', 'mask', 'q'],
                                output_names=['form_str'], function=fdr), 
                       name='form_neg', iterfield='p_im', nested=True)
    form_neg.inputs.p_im = p_file_neg
    form_neg.inputs.mask = mask
    form_neg.inputs.q = cor['q']
    neg = form_neg.run().outputs.form_str
    
    corrected_pos = MapNode(ImageMaths(suffix='_fdr'), name='corrected_pos', iterfield=['in_file', 'op_string'], nested=True)
    corrected_pos.inputs.in_file = p_file_pos
    corrected_pos.inputs.op_string = pos
    corrected_pos = corrected_pos.run().outputs.out_file
    
    corrected_neg = MapNode(ImageMaths(suffix='_fdr'), name='corrected_neg', iterfield=['in_file', 'op_string'], nested=True)
    corrected_neg.inputs.in_file = p_file_neg
    corrected_neg.inputs.op_string = neg
    corrected_neg = corrected_neg.run().outputs.out_file
    
    flip = MapNode(ImageMaths(op_string='-mul -1'), name='corrected_neg', iterfield=['in_file'], nested=True)
    flip.inputs.in_file = corrected_neg
    corrected_neg = flip.run().outputs.out_file
    
    corrected_all = MapNode(BinaryMaths(operation='add'), name='corrected_all', iterfield=['in_file', 'operand_file', 'out_file'], nested=True)
    corrected_all.inputs.in_file = corrected_pos
    corrected_all.inputs.operand_file = corrected_neg
    
    out_names = []
    for i in range(len(corrected_pos)):
        out_names.append(os.getcwd()+'/fdr_corrected_{i}.nii.gz'.format(i=i))
        
    corrected_all.inputs.out_file = out_names
    
    return corrected_all.run().outputs.out_file

def FWE(zstat, mask, cor):
    from updated.correction.functions import fwe, neg
    from nipype import MapNode, Function
    from nipype.interfaces.fsl import SmoothEstimate, Threshold, BinaryMaths
    import os
    
    smoothness = MapNode(SmoothEstimate(mask_file=mask), 
                         name='smoothness', iterfield=['zstat_file'], nested=True)
    smoothness.inputs.zstat_file = zstat
    smoothness = smoothness.run().outputs.resels
    
    thresh = MapNode(Function(input_names=['p', 'resels'],
                              output_names=['thresh'], function=fwe), 
                     name='thresh', iterfield=['resels'], nested=True)
    thresh.inputs.p = cor['p']
    thresh.inputs.resels = smoothness
    thresh = thresh.run().outputs.thresh
    
    fwe_nonsig0 = MapNode(Threshold(direction='above'), name='fwe_nonsig0', iterfield=['in_file', 'thresh'], nested=True)
    fwe_nonsig0.inputs.in_file = zstat
    fwe_nonsig0.inputs.thresh = thresh
    pos = fwe_nonsig0.run().outputs.out_file
    
    fwe_nonsig1 = MapNode(Threshold(direction='below'), name='fwe_nonsig1', iterfield=['in_file', 'thresh'], nested=True)
    fwe_nonsig1.inputs.in_file = pos
    fwe_nonsig1.inputs.thresh = neg(thresh)
    neg = fwe_nonsig1.run().outputs.out_file
    
    fwe_thresh = MapNode(BinaryMaths(operation='sub'), name='fwe_thresh', iterfield=['in_file', 'operand_file', 'out_file'], nested=True)
    fwe_thresh.inputs.in_file = zstat
    fwe_thresh.inputs.operand_file = neg
    
    out_names = []
    for i in range(len(zstat)):
        out_names.append(os.getcwd()+'/fwe_corrected_{i}.nii.gz'.format(i=i))
        
    fwe_thresh.inputs.out_file = out_names
    
    return fwe_thresh.run().outputs.out_file

def clusterFWE(zstat, copes, mask, cor):
    #RETURNS CORRECTED 1-P IMAGE
    from nipype import MapNode
    from nipype.interfaces.fsl import SmoothEstimate, Cluster, BinaryMaths
    import os
    
    smoothness = MapNode(SmoothEstimate(mask_file=mask), 
                         name='smoothness', iterfield=['zstat_file'], nested=True)
    smoothness.inputs.zstat_file = zstat
    smoothness = smoothness.run().outputs
    
    cluster_pos = MapNode(Cluster(out_index_file=True, out_localmax_txt_file=True, out_threshold_file=True), 
                          name='cluster_pos', iterfield=['in_file', 'cope_file', 'volume', 'dlh'], nested=True)
    cluster_pos.inputs.connectivity = cor['connectivity']
    cluster_pos.inputs.threshold = cor['zthreshold']
    cluster_pos.inputs.pthreshold = cor['pthreshold']
    cluster_pos.inputs.volume = smoothness.volume
    cluster_pos.inputs.dlh = smoothness.dlh
    cluster_pos.inputs.in_file = zstat
    cluster_pos.inputs.cope_file = copes
    cluster_pos = cluster_pos.run().outputs.threshold_file
    
    zstat_inv = MapNode(BinaryMaths(operation='mul', operand_value=-1), 
                        name='zstat_inv', iterfield=['in_file'], nested=True)
    zstat_inv.inputs.in_file = zstat
    zstat_inv = zstat_inv.run().outputs.out_file
    
    cluster_neg = MapNode(Cluster(out_index_file=True, out_localmax_txt_file=True, out_threshold_file=True), 
                          name='cluster_neg', iterfield=['in_file', 'cope_file', 'volume', 'dlh'], nested=True)
    cluster_neg.inputs.connectivity = cor['connectivity']
    cluster_neg.inputs.threshold = cor['zthreshold']
    cluster_neg.inputs.pthreshold = cor['pthreshold']
    cluster_neg.inputs.volume = smoothness.volume
    cluster_neg.inputs.dlh = smoothness.dlh
    cluster_neg.inputs.in_file = zstat_inv
    cluster_neg.inputs.cope_file = copes
    cluster_neg = cluster_neg.run().outputs.threshold_file
    
    
    cluster_inv = MapNode(BinaryMaths(operation='mul', operand_value=-1), 
                          name='cluster_inv', iterfield=['in_file'], nested=True)
    cluster_inv.inputs.in_file = cluster_neg
    cluster_inv = cluster_inv.run().outputs.out_file
        
    cluster_all = MapNode(BinaryMaths(operation='add'), 
                          name='cluster_all', iterfield=['in_file', 'operand_file', 'out_file'], nested=True)
    cluster_all.inputs.in_file = cluster_pos
    cluster_all.inputs.operand_file = cluster_inv
    
    out_names = []
    for i in range(len(cluster_pos)):
        out_names.append(os.getcwd()+'/cluster_corrected_{i}.nii.gz'.format(i=i))
        
    cluster_all.inputs.out_file = out_names
    
    cluster_all = cluster_all.run().outputs.out_file
    
    return cluster_all


