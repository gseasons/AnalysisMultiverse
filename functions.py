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

def parse_xml(xml, goal, mask):
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
    
    def get_warps(self):
        from nipype.interfaces.fsl import InvWarp
        genwarps = Workflow('genwarps')
        genwarps.base_dir = os.getcwd() #self.base_dir
        
        inputnode = Node(IdentityInterface(fields=['brain', 'ref_file']), name='inputnode')
        outnode = Node(IdentityInterface(fields=['warp', 'invwarp']), name='outnode')
        
        prelim = Node(FLIRT(dof=12, output_type='NIFTI_GZ'), name='prelim')
        warp = Node(FNIRT(field_file=True), name='warp')
        
        warp_inv = Node(InvWarp(), name='warp_inv')
        
        genwarps.connect([(inputnode, prelim, [('brain', 'in_file')]),
                          (inputnode, prelim, [('ref_file', 'reference')]),
                          (inputnode, warp, [('brain', 'in_file')]),
                          (inputnode, warp, [('ref_file', 'ref_file')]),
                          (prelim, warp, [('out_matrix_file', 'affine_file')]),
                          (warp, warp_inv, [('field_file', 'warp')]),
                          (inputnode, warp_inv, [('brain', 'reference')]),
                          (warp, outnode, [('field_file', 'warp')]),
                          (warp_inv, outnode, [('inverse_warp', 'invwarp')]),
                          ])
        
        return genwarps
        
    def apply_warps(self):
        appwarps = Workflow('appwarps')
        appwarps.base_dir = os.getcwd()
        
        inputnode = Node(IdentityInterface(fields=['feat_dir', 'warp_file', 'ref_file', 'needwarp']), name='inputnode')
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
        
        applywarp_c = MapNode(ApplyWarp(), name='applywarp_c', iterfield=['in_file'])
        applywarp_v = MapNode(ApplyWarp(), name='applywarp_v', iterfield=['in_file'])
        applywarp_bold = Node(ApplyWarp(), name='applywarp_bold')
        
        appwarps.connect([(inputnode, applywarp_c, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_v, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_bold, [('ref_file', 'ref_file')]),
                          (inputnode, applywarp_c, [('warp_file', 'field_file')]),
                          (inputnode, applywarp_v, [('warp_file', 'field_file')]),
                          (inputnode, applywarp_bold, [('warp_file', 'field_file')]),
                          (inputnode, selectfiles, [('feat_dir', 'base_directory')]),
                          (inputnode, ret, [('needwarp', 'needwarp')]),
                          (selectfiles, ident, [('cope', 'cope')]),
                          (selectfiles, ident, [('varcope', 'varcope')]),
                          (selectfiles, ident, [('bold', 'bold')]),
                          (inputnode, ident, [('needwarp', 'needwarp')]),
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
                          (ret, outnode, [('bold', 'bold')]),
                          ])
        
        return appwarps