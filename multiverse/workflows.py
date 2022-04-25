#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:35:10 2021

@author: grahamseasons
"""
def write_out(base_dir, pipeline_st, task):
    outputs = list(locals().keys())[3:]
    from nipype.interfaces import DataSink
    from nipype import Node
    import re, os, glob
    sink = DataSink(base_directory=base_dir, parameterization=True)
    folder_name = task
    test_str = vars()[outputs[0]]
    if isinstance(test_str, list):
        test_str = test_str[0]
        if isinstance(test_str, list):
            test_str = test_str[0]
    session = re.search('/(_sessions_[A-Za-z0-9]+)/', test_str)
    run = re.search('/(_runs_[0-9]+)/', test_str)
    #subject = re.search('/(_subject_[0-9A-Za-z]+/)', test_str).group(1)
    
    if session:
        folder_name += '/' + session.group(1) + '/'
    if run:
        folder_name += '/' + run.group(1) + '/'
    if folder_name == task:
        folder_name += '/'
    
    #folder_name += wf_name + '/' + subject
    
    for i, out in enumerate(outputs):
        if out:
            
            setattr(sink.inputs, 'pipelines/' + task + '.@' + str(i), vars()[out])
            
    sink.run()