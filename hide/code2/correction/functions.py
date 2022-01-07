#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 12:46:46 2021

@author: grahamseasons
"""
from functions import insert
from nipype.utils.functions import getsource
from workflows import write_out
import re

def fdr(p_im, mask, q):
    from nipype.interfaces.base import CommandLine
    import re
    cmd = ('fdr -i {p_im} -m {mask} -q {q}')
    cl = CommandLine(cmd.format(p_im=p_im, mask=mask, q=q))
    results = cl.run().runtime.stdout
    thresh = 1 - float(re.search('.*([0-9\.]+)', results).group(0))
    
    form = '-mul -1 -add 1 -thr {thresh}'.format(thresh=thresh)
    
    return form

def fwe(p, resels):
    from nipype.interfaces.base import CommandLine
    import re
    p_two = p / 2
    cmd = ('ptoz {p} -2 -g {resels}')
    cl = CommandLine(cmd.format(p=p_two, resels=resels))
    results = cl.run().runtime.stdout
    thresh = float(re.search('([0-9\.]+)', results).group(1))
    
    return thresh

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

def correction(zstat, copes, mask, cor):
    from correction.workflows import FDR, FWE, clusterFWE
    method = cor['method']
    if method == 'fdr':
        corrected = FDR(zstat, mask, cor)
    elif method == 'fwe':
        corrected = FWE(zstat, mask, cor)
    elif method == 'clust':
        corrected = clusterFWE(zstat, copes, mask, cor)
    else:
        print('ERROR: CORRECTION METHOD NOT IMPLEMENTED')
        
    return corrected

def get_sink(inputs):
    func_str = getsource(write_out)
    ind = func_str.find('):')
    params = ', ' + ', '.join(inputs)
    func_str = insert(func_str, ind, params)
        
    return func_str
        
        