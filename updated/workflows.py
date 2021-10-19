#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:35:10 2021

@author: grahamseasons
"""
def write_out(base_dir, pipeline_st, task):
    outputs = list(locals().keys())[5:]
    from nipype.interfaces import DataSink
    from nipype import Node
    import re, os
    sink = Node(DataSink(base_directory=base_dir, parameterization=True), name='sink')
    folder_name = task
    session = re.search('/(_sessions_[A-Za-z0-9]+)/', vars()[outputs[0]])
    run = re.search('/(_runs_[0-9]+)/', vars()[outputs[0]])
    subject = re.search('/(_subject_[0-9A-Za-z]+/)', vars()[outputs[0]]).group(1)
    
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