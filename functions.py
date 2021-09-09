#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 11:17:17 2021

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

from nipype.interfaces.base import CommandLine
from nipype.interfaces.fsl import Cluster, BinaryMaths, Threshold, SmoothEstimate
import re

def get_wm(files):
    return files[-1]

def get_bright_thresh(medianval):
    return 0.75 * medianval

def getthreshop(thresh):
    return '-thr %.10f -Tmin -bin' % (0.1 * thresh[1])

def FWE(base_dir, zstat, p, mask):
        #RETURNS CORRECTED ZSCORE IMAGE
        FWE = Workflow('FWE')
        FWE.base_dir = base_dir
        
        inputnode = Node(IdentityInterface(fields=['zstat', 'mask', 'p']), name='inputnode')
        inputnode.inputs.zstat = zstat
        inputnode.inputs.p = p
        inputnode.inputs.mask = mask
        outnode = Node(IdentityInterface(fields=['corrected']), name='outnode')
        
        smoothness = MapNode(SmoothEstimate(), name='smoothness', iterfield=['zstat_file'], nested=True)
        
        def fwe(p, resels):
            from nipype.interfaces.base import CommandLine
            import re
            p_two = p / 2
            cmd = ('ptoz {p} -2 -g {resels}')
            cl = CommandLine(cmd.format(p=p_two, resels=resels))
            results = cl.run().runtime.stdout
            thresh = float(re.search('([0-9\.]+)', results).group(1))
            
            return thresh
        
        thresh = MapNode(Function(input_names=['p', 'resels'],
                               output_names=['thresh'], function=fwe), name='thresh', iterfield=['resels'], nested=True)
        
        def neg(vals):
            isnested = any(isinstance(i, list) for i in vals)
            if isnested:
                val = [[-1 * val for val in sublst] for sublst in vals]
            else:
                if type(vals) == list:
                    val = [-1 * val for val in vals]
                else:
                    val = -1 * val
                
            return val
        
        fwe_nonsig0 = MapNode(Threshold(direction='above'), name='fwe_nonsig0', iterfield=['in_file', 'thresh'], nested=True)
        fwe_nonsig1 = MapNode(Threshold(direction='below'), name='fwe_nonsig1', iterfield=['in_file', 'thresh'], nested=True)
        fwe_thresh = MapNode(BinaryMaths(operation='sub'), name='fwe_thresh', iterfield=['in_file', 'operand_file'], nested=True)
        
        FWE.connect([(inputnode, smoothness, [('zstat', 'zstat_file'),
                                              ('mask', 'mask_file')]),
                     (inputnode, thresh, [('p', 'p')]),
                     (inputnode, fwe_thresh, [('zstat', 'in_file')]),
                     (inputnode, fwe_nonsig0, [('zstat', 'in_file')]),
                     (fwe_nonsig0, fwe_nonsig1, [('out_file', 'in_file')]),
                     (smoothness, thresh, [('resels', 'resels')]),
                     (thresh, fwe_nonsig0, [('thresh', 'thresh')]),
                     (thresh, fwe_nonsig1, [(('thresh', neg), 'thresh')]),
                     (fwe_nonsig1, fwe_thresh, [('out_file', 'operand_file')]),
                     (fwe_thresh, outnode, [('out_file', 'corrected')]),
                     ])
        
        return FWE
        
        
        #TODO: NON CLUSTER FWE, GENERATE ANALYSIS WORKFLOW, RENAME OUTPUTS TO USER FRIENDLY NAMES
        
def clusterFWE(base_dir, p, z_thresh, mask, connectivity, copes, zstat):
        #RETURNS CORRECTED 1-P IMAGE
        from nipype import IdentityInterface
        cluster = Workflow('Cluster')
        cluster.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['zstat', 'mask', 'connectivity', 'copes' 
                                                   'threshold', 'pthreshold']), name='inputnode')
        inputnode.inputs.pthreshold = p
        inputnode.inputs.threshold = z_thresh
        inputnode.inputs.mask = mask
        inputnode.inputs.connectivity = connectivity
        inputnode.inputs.copes = copes
        inputnode.inputs.zstat = zstat
        
        outnode = Node(IdentityInterface(fields=['corrected']), name='outnode')
        
        smoothness = MapNode(SmoothEstimate(), name='smoothness', iterfield=['zstat_file'], nested=True)
        
        cluster_pos = MapNode(Cluster(out_index_file=True, out_localmax_txt_file=True, out_threshold_file=True), 
                              name='cluster_pos', iterfield=['in_file', 'cope_file', 'volume', 'dlh'], nested=True)
        cluster_neg = MapNode(Cluster(out_index_file=True, out_localmax_txt_file=True, out_threshold_file=True), 
                              name='cluster_neg', iterfield=['in_file', 'cope_file', 'volume', 'dlh'], nested=True)
        
        zstat_inv = MapNode(BinaryMaths(operation='mul', operand_value=-1), name='zstat_inv', iterfield=['in_file'], nested=True)
        cluster_inv = MapNode(BinaryMaths(operation='mul', operand_value=-1), name='cluster_inv', iterfield=['in_file'], nested=True)
        cluster_all = MapNode(BinaryMaths(operation='add'), name='cluster_all', iterfield=['in_file', 'operand_file'], nested=True)
        
        cluster.connect([(inputnode, smoothness, [('zstat', 'zstat_file'),
                                                  ('mask', 'mask_file')]),
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
        
        return cluster
        
#TODO: FIX THIS
def FDR(base_dir, zstat, p, mask):
        #RETURNS CORRECTED 1-P IMAGE
        #SEPARATE POSITIVE AND NEGATIVE FORMS
        FDR = Workflow('FDR')
        FDR.base_dir = base_dir #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['zstat', 'mask', 'p']), name='inputnode')
        inputnode.inputs.zstat = zstat
        inputnode.inputs.p = p
        inputnode.inputs.mask = mask
        outnode = Node(IdentityInterface(fields=['corrected']), name='outnode')
        
        p_file = MapNode(ImageMaths(op_string='-ztop', suffix='_pval'), name='p_file', iterfield=['in_file'], nested=True)
        
        def fdr(p_im, mask, q): #q = 0.05
            from nipype.interfaces.base import CommandLine
            import re
            cmd = ('fdr -i {p_im} -m {mask} -q {q}')
            cl = CommandLine(cmd.format(p_im=p_im, mask=mask, q=q))
            results = cl.run().runtime.stdout
            thresh = 1 - float(re.search('([0-9\.]+)', results).group(1))
            
            form = '-mul -1 -add 1 -thr {thresh} -mas {mask}'.format(thresh=thresh, mask=mask)
            
            return form
        
        form = MapNode(Function(input_names=['p_im', 'mask', 'q'],
                               output_names=['form_str'], function=fdr), name='form', iterfield=['p_im'], nested=True)
        
        
        corrected = MapNode(ImageMaths(suffix='_fdr'), name='corrected', iterfield=['in_file', 'op_string'], nested=True)
        
        FDR.connect([(inputnode, p_file, [('zstat', 'in_file')]),
                     (inputnode, form, [('mask', 'mask'),
                                        ('p', 'q')]),
                     (p_file, form, [('out_file', 'p_im')]),
                     (form, corrected, [('form_str', 'op_string')]),
                     (p_file, corrected, [('out_file', 'in_file')]),
                     (corrected, outnode, [('out_file', 'corrected')])
                     ])
        
        return FDR