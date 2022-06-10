#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 10:28:41 2021

@author: grahamseasons
"""
from nipype.utils.functions import getsource
from functions import insert
from workflows import write_out
import re

def get_bright_thresh(medianval):
    """Calculate brightness threshold for SUSAN workflow"""
    return 0.75 * medianval

def getthreshop(thresh):
    """Generate mask for SUSAN workflow"""
    return '-thr %.10f -Tmin -bin' % (0.1 * thresh[1])

def get_wm(files):
    """Return white matter file from FAST"""
    return files[-1]

def decision(mask, mc_mean, mc, st, slice_correct='', mean_vol=''):
    """Ability to turn on/off slice timing correction as well as select median volume if not corrected to mean"""
    from nipype.interfaces.fsl import Threshold, ExtractROI, FLIRT, ImageStats, UnaryMaths
    import nibabel as nib
    import os, re
    from shutil import copy2
    resample = False
    if nib.load(mc).shape[0:-1] != nib.load(mask).shape:
        resample = True
    
    if 'mean_reg' in mean_vol and mean_vol.count('.nii.gz') == 2:
        new_name = re.search('(sub-[0-9A-Za-z\-_]+)', mean_vol).group(1) + '_mean_reg.nii.gz'
        mean_vol = copy2(mean_vol, os.path.join(os.getcwd(), new_name))
        
    mask = Threshold(in_file=mask, thresh=0, args='-bin').run().outputs.out_file
    mask = UnaryMaths(in_file=mask, operation='fillh').run().outputs.out_file
    
    if mc_mean and slice_correct:
        if resample:
            mask = FLIRT(in_file=mask, reference=mean_vol, apply_xfm=True, uses_qform=True, interp='nearestneighbour').run().outputs.out_file
            if not int(ImageStats(in_file=mask, op_string='-R').run().outputs.out_stat[1]):
                mask = ''
        return mean_vol, st, mask
    elif mc_mean:
        if resample:
            mask = FLIRT(in_file=mask, reference=mean_vol, apply_xfm=True, uses_qform=True, interp='nearestneighbour').run().outputs.out_file
            if not int(ImageStats(in_file=mask, op_string='-R').run().outputs.out_stat[1]):
                mask = ''
        return mean_vol, mc, mask
    elif slice_correct:
        size = nib.load(mc).shape
        get_mid = ExtractROI(in_file=mc, t_min=round(size[-1]/2), t_size=1)
        mc_m = get_mid.run().outputs.roi_file
        if resample:
            mask = FLIRT(in_file=mask, reference=mc_m, apply_xfm=True, uses_qform=True, interp='nearestneighbour').run().outputs.out_file
            if not int(ImageStats(in_file=mask, op_string='-R').run().outputs.out_stat[1]):
                mask = ''
        return mc_m, st, mask
    else:
        size = nib.load(mc).shape
        get_mid = ExtractROI(in_file=mc, t_min=round(size[-1]/2), t_size=1)
        mc_m = get_mid.run().outputs.roi_file
        if resample:
            mask = FLIRT(in_file=mask, reference=mc_m, apply_xfm=True, uses_qform=True, interp='nearestneighbour').run().outputs.out_file
            if not int(ImageStats(in_file=mask, op_string='-R').run().outputs.out_stat[1]):
                mask = ''
        return mc_m, mc, mask
    
    
def strip_container(in_file):
    """Remove container"""
    if type(in_file) == list:
        return in_file[0]
    else:
        return in_file

def function_str(name, dic=''):
    """Injects code which allows parameters to be dynamically assigned to workflows"""
    from preprocessing.workflows import registration, smooth, regress, mni
    valid_functions = ['registration', 'smooth', 'regress', 'mni']
    if name in valid_functions:
        func_str = getsource(vars()[name])
        try: 
            out = []
            for names in dic[name].keys():
                out.append(re.search('[A-Za-z]+_([A-Za-z_]+)', names).group(1))
                
            ind = func_str.find('):')
            params = ', ' + ', '.join(out)
            func_str = insert(func_str, ind, params)
            out = [element for element in out if '_' in element]
            
            workflowfind = list(re.finditer('\n(\n)(\s+)[A-Za-z]+.run()', func_str))
            if not workflowfind:
                workflowfind = list(re.finditer('\n(\n)(\s+)[A-Za-z\s=]+.run()', func_str))
            
            for search in reversed(workflowfind):
                ind = search.start(1)
                block = '\n' + search.group(2) + (search.group(2) + search.group(2)).join(["for param in {params}:\n",
                "search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)\n",
                "if vars()[param]: setattr(vars()[search.group(1)].inputs, search.group(2), vars()[param])\n",
                "else: setattr(vars()[search.group(1)].inputs, search.group(2), Undefined)\n"])
                
                func_str = insert(func_str, ind, block.format(params=out))
            return func_str, re.search('def ' + name + '\(([A-Za-z_,0-9\s]+)\)', func_str).group(1).split(', ')
        except:
            return func_str, re.search('def ' + name + '\(([A-Za-z_,0-9\s]+)\)', func_str).group(1).split(', ')
        
def get_sink(inputs):
    """Expands datasink to include inputs"""
    func_str = getsource(write_out)
    ind = func_str.find('):')
    params = ', ' + ', '.join(inputs)
    func_str = insert(func_str, ind, params)
    return func_str
