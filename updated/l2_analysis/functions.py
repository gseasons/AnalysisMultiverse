#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 11:45:06 2021

@author: grahamseasons
"""
from updated.functions import insert
from nipype.utils.functions import getsource
from updated.workflows import write_out
import re

def groupscans(copes, varcopes):
    outcopes = []
    outvarcopes = []
    multiple = type(copes[0])
    if multiple == list:
        numcon = len(copes[0])
    else:
        return [copes], [varcopes], len(copes)
        
    for i in range(numcon):
        coperuns = [run[i] for run in copes]
        varcoperuns = [run[i] for run in varcopes]
        outcopes.append(coperuns)
        outvarcopes.append(varcoperuns)
        
    return outcopes, outvarcopes, numcon

def get_sink(inputs, files):
    if type(files) != list:
        files = [files]
        
    func_str = getsource(write_out)
    ind = func_str.find('):')
    params = ', ' + ', '.join(inputs)
    func_str = insert(func_str, ind, params)
    
    search = re.search('\n(\n)(\s+)(setattr)', func_str)
    ind = search.start(1)
    
    block = '\n' + search.group(2) + (search.group(2) + search.group(2)).join(["if isinstance(vars()[out], str) and os.path.isdir(vars()[out]):\n",
        "for file in {files}:\n",
        "    file_out = vars()[out] + '/' + file\n",
        "    if os.path.isdir(file_out) or os.path.isfile(file_out):\n",
        "        setattr(sink.inputs, 'pipelines/' + task + '.@' + str(i) + file, file_out)\n"])
        
    func_str = insert(func_str, search.start(3), "\n    "+search.group(2))
    func_str = insert(func_str, search.start(3), "else:")
    func_str = insert(func_str, ind, block)
        
    return func_str.format(files=str(files))