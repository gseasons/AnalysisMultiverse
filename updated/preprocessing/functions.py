#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 10:28:41 2021

@author: grahamseasons
"""
from nipype.utils.functions import getsource
import re

def get_bright_thresh(medianval):
    return 0.75 * medianval

def getthreshop(thresh):
    return '-thr %.10f -Tmin -bin' % (0.1 * thresh[1])

def insert(string, ind, new):
    return string[:ind] + new + string[ind:]

def get_wm(files):
    return files[-1]

def decision(mc_mean, mc, st, slice_correct, mean_vol):
    if mean_vol and slice_correct:
        return mc_mean, st
    elif mean_vol:
        return mc_mean, mc
    elif slice_correct:
        return mc, st
    else:
        return mc, mc

def function_str(name, dic=''):   
    from updated.preprocessing.workflows import registration, smooth
    valid_functions = ['registration', 'smooth']
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
                
            for search in reversed(list(re.finditer('\n(\n)(\s+)[A-Za-z]+.run()', func_str))):
                ind = search.start(1)
                block = '\n' + search.group(2) + (search.group(2) + search.group(2)).join(["for param in {params}:\n",
                "search = re.search('([A-Za-z]+)_([A-Za-z_]+)', param)\n",
                "setattr(vars[search.group(1)].inputs, search.group(2), vars()[param])\n"])
                
                func_str = insert(func_str, ind, block).format(params=out)
            return func_str, re.search('def ' + name + '\(([A-Za-z_,0-9\s]+)\)', func_str).group(1).split(', ')
        except:
            return func_str, re.search('def ' + name + '\(([A-Za-z_,0-9\s]+)\)', func_str).group(1).split(', ')
        
def out(base_dir, pipeline_st, task, param_dic, slicetiming, smooth, mc_par, mc_img, mc_mean, brain, outliers, plots, out_mat, warped, segment, warp_file):
    outputs = list(locals().keys())[3:]
    from nipype.interfaces import DataSink
    from nipype import Node
    import re, os
    sink = Node(DataSink(base_directory=base_dir, parameterization=False), name='sink')
    folder_name = task
    session = re.search('/(_sessions_[A-Za-z0-9]+)/', vars()[outputs[0]])
    run = re.search('/(_runs_[0-9]+)/', vars()[outputs[0]])
    subject = re.search('/(_subject_[0-9A-Za-z]+/)', vars()[outputs[0]]).group(1)
    pipeline = re.search('/mapflow/[A-Za-z_]+([0-9]+)', os.getcwd()).group(1)
    pipeline = 'pipeline_' + str(int(pipeline) + pipeline_st)
    
    sink.inputs.container = pipeline
    
    if session:
        folder_name += '/' + session.group(1) + '/'
    if run:
        folder_name += '/' + run.group(1) + '/'
    if folder_name == task:
        folder_name += '/'
            
    folder_name += 'preprocess/' + subject
    
    outputs = []
    
    for i, out in enumerate(outputs):
        setattr(sink.inputs, folder_name + '.@' + str(i), vars()[out])
            
    sink.run()
