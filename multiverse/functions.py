#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:54:23 2021

@author: grahamseasons
"""
from nipype import Node, IdentityInterface, Function
from nipype.interfaces.base import Undefined
import re, random
from collections import Counter
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
import os

exp_dir = '/scratch'
out_dir = exp_dir + '/processed'

def no_mask(file):
    raise FileNotFoundError("Specified mask '{m}' does not exist".format(m=file))
    
def invalid(mapped):
    raise SyntaxError("Input paramater '{name}' in an unsupported format. Acceptable formats are (brain region, thr), (brain region, thr_min, thr_max, or path_to_mask".format(name=mapped))

def generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle, ind_start, frame):
    """Generates dictionaries of analysis parameters for workflow from GA output.
       i.e. Converts numeric values into values accepted by nodes, as defined in links"""
    container = np.zeros((1, pop.shape[0]), str)
    indexes = [(ind_start, ind_start+pop.shape[0])]
    #Many of the below are accessed through var() -> they are in use, do not delete
    expand_inputs = {}
    
    preprocess = {}
    level1 = {}
    
    unique_cols = []
    unique_cols_temp = []
    unique_vals = []
    
    if multiscan:
        level2 = {}
        
    level3 = {}
    correction = {}
    dic = 'preprocess'
    valid = []
    master = {}
    previous = ''
    counter = 0
    for key in map_genes:
        gene = map_genes[key]
        keys = list(gene.keys())
        node_name = re.search('([A-Za-z0-9]+)_', keys[0]).group(1)
        if not previous or 'end' in previous:
            if key - 1 in map_genes:
                get_value = key - 1
            else:
                get_value = key
            #Checks to see if multiple scans per subject, if not don't include second level in analysis
            if not multiscan and list(map_genes[get_value].values())[0] == 'level2':
                lock = True
                counter += 1
                continue
            elif vars().get('lock', False):
                if node_name == 'end':
                    lock = False
                    dic = gene[keys[0]]
                    continue
                else:
                    continue
            previous = node_name
            
        if previous != node_name:
            #Check to see if any nodes need to be added in from links json
            for link in links[dic]:
                if previous in link[:len(previous)]:
                    connect = links[dic][link]
                    #Check to see if direct copy, or formatted copy of node
                    if type(connect) == list:
                        group = re.search('([A-Za-z0-9]+)_', connect[0]).group(1)
                        if group not in vars()[dic]:
                            for opt in valid:
                                if group in vars()[opt]:
                                    break
                        vals = vars()[dic][group][connect[0]]
                        check = vars()[dic][group][connect[1]]
                        #Assign mutually exclusive values
                        if len(connect) == 3:
                            rule = connect[2]
                            vars()[dic][previous][link] = [val if check[c] != rule else Undefined for c, val in enumerate(vals)]
                        else:
                            vars()[dic][previous][link] = [val if check[c] else Undefined for c, val in enumerate(vals)]
                    else:
                        #Direct copy
                        try:
                            group = re.search('([A-Za-z0-9]+)_', connect).group(1)
                            if group in vars()[dic]:
                                vars()[dic][previous][link] = vars()[dic][group][connect]
                            else:
                                for opt in valid:
                                    if group in vars()[opt]:
                                        break
                                vars()[dic][previous][link] = vars()[opt][group][connect]
                        except:
                            raise SyntaxError("Node name {pre} in links violates naming convention. Please keep name to alphanumeric characters.".format(pre=connect))
            #Dynamic input dictionary   
            if 'F' == previous[0]:
                expand_inputs[previous[1:]] = vars()[dic][previous]
                
            previous = node_name
        #Check if analysis stage complete, define where pipelines split from each other
        if 'end' in keys[0]:
            container, pholder, indexes = define_paths(container, vars()[dic], indexes)
            vars()[dic+'_old'] = vars()[dic].copy()
            valid.append(dic+'_old')
            vars()[dic].clear()
            vars()[dic].update(pholder)
            master[dic] = vars()[dic]
            dic = gene[keys[0]]
            counter += 1
            continue
        #Create analysis stage dictionary
        if node_name not in vars()[dic]:
            vars()[dic][node_name] = {}
            
        values = params[key-counter,:]
        #Convert to proper datatypes (float to int that can be int)
        isint = False
        if values.dtype == float:
            check_vals = [val.is_integer() for val in values]
            if sum(check_vals) == len(check_vals):
                isint = True
        #Skip linked nodes (dealt with above)
        if keys[0][0] == '!':
            continue
        #non-numeric parameters mapped, or constructing a dictionary of parameters
        if len(gene) > 1 or '~construct~' in keys[0]:
            for l, i in enumerate(values):
                if round(i) in gene:
                    mapped = gene[round(i)]
                else:
                    mapped = i
                #special case of a dictionary containing mutually exclusive parameters (see applywarp_ XFM Matrix)
                #expands into multiple parameters
                if type(mapped) == dict and keys[0][-1] == '_':
                    for k in mapped:
                        if keys[0][-1] == '_':
                            param = keys[0] + k
                        else:
                            param = node_name + '_' + k
                        if param not in vars()[dic][node_name]:
                             vars()[dic][node_name][param] = []
                        #Add parameter name to list of independently chosen parameters
                        if param not in unique_cols and param not in unique_cols_temp:
                            unique_cols_temp.append(param)
                        vars()[dic][node_name][param].append(mapped[k])
                #construct dictionary
                elif '~construct~' in keys[0]:
                    var_name = re.search('_([A-Za-z]+)', keys[0]).group(1)
                    key_name = re.search('_([A-Za-z]+)$', keys[0]).group(1)
                    param = node_name + '_' + var_name
                    
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    #check for every value being an integer
                    if not isint:
                        rand = random.Random(i)
                        #list of networks targeted with atlas (form of a tuple, if 1 threshold given creates range of possible thresholds, if two uses as range)
                        if type(mapped) == list:
                            mapped = [(m[0], rand.randint(int(m[1])-wiggle if int(m[1])-wiggle > 0 else 0, int(m[1])+wiggle if int(m[1])+wiggle > 95 else 95)) if len(m) == 2 and isinstance(m, (tuple,list)) else (m[0], rand.randint(int(m[1]), int(m[2]))) if isinstance(m, (tuple,list)) else m if isinstance(m, str) else invalid(param+'_'+key_name) for m in mapped]
                        #tuple of integers generate from range
                        elif mapped == tuple: #as far as i can tell this is redundant/never used
                            mapped = rand.randint(mapped[0], mapped[1])
                        else:
                            mapped = i
                    else:
                        if isinstance(mapped, float):
                            mapped = int(mapped)
                    #check to see if all values mapped (list of dictionaries to construct -> verifies dictionary added for each pipeline)      
                    if len(vars()[dic][node_name][param]) != pop.shape[0]:
                        construction_key = key_name
                        vars()[dic][node_name][param].append({key_name: mapped})
                    else:
                        #check to make sure parameter is supposed to be included in dictionary (i.e. ROI coords only relevant for ROI seed based analysis, not atlas or data driven)
                        if key_name in links[dic]:
                            if vars()[dic][node_name][param][l][construction_key] not in links[dic][key_name]:
                                continue
                        #add parameter to correct slot in list
                        vars()[dic][node_name][param][l][key_name] = mapped
                else:
                    #add mapped values (no special cases)
                    param = keys[0]
                        
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    
                    vars()[dic][node_name][param].append(mapped)

            if param not in unique_cols and param not in unique_cols_temp:
                unique_cols_temp.append(param)
        else:
            #convert float values that are supposed to be integers into integers
            if values.dtype == float:
                if isint:
                    values = values.astype(int)
            #add values to dictionary
            vars()[dic][node_name][keys[0]] = values
            if keys[0] not in unique_cols and keys[0] not in unique_cols_temp:
                unique_cols_temp.append(keys[0])
        
        #get unique values in list to later construct dataframe
        for param in unique_cols_temp:
            node_name_temp = re.search('([A-Za-z0-9]+)_', param).group(1)
            try:
                unique_vals.append(vars()[dic][node_name_temp][param])
            except KeyError:
                for opt in valid:
                    if node_name_temp in vars()[opt]:
                        break
                unique_vals.append(vars()[opt][node_name_temp][param])
                
        unique_cols += unique_cols_temp
        unique_cols_temp = []
        
    #create dataframe of pipelines
    unique_cols += ['R', 'P', 'Score']
    l = len(unique_vals[-1])
    unique_vals += [[0]*l, [0]*l, [0]*l]
    
    pipeline_def = pd.DataFrame(data=np.array(unique_vals, dtype=object).transpose(), columns=unique_cols)
    #if dataframe already exists, add to it (for genetic algorithm)
    if type(frame) != str:
        pipeline_def = pd.concat([frame, pipeline_def], ignore_index=True)

    return master, expand_inputs, pipeline_def
    
def metadata(filename, data):
    """Grab TR time from metadata"""
    from bids.layout import BIDSLayout
    layout = BIDSLayout(data)
    meta_data = layout.get_metadata(filename)
    TR = meta_data['RepetitionTime']
    return TR

def event_grabber(file, data):
    """Get event file/information for subject"""
    import re
    from bids.layout import BIDSLayout
    layout = BIDSLayout(data)
    
    task = re.search('task-([0-9A-Za-z]+)_bold', file).group(1)
    
    if 'rest' in task:
        return ''
    
    event_file = layout.get(task=task, extension='.tsv')
    
    if len(event_file) > 1:
        sub = re.search('/sub-([0-9A-Za-z]+)/', file).group(1)
        ses = re.search('_ses-([A-Za-z]+)_task', file)
        run = re.search('run-([0-9]+)', file)
        
        if ses and run:
            event_file = layout.get(task=task, session=ses.group(1), run=run.group(1), subject=sub, extension='.tsv')
        elif ses:
            event_file = layout.get(task=task, session=ses.group(1), subject=sub, extension='.tsv')
        elif run:
            event_file = layout.get(task=task, run=run.group(1), subject=sub, extension='.tsv')
    elif not len(event_file):
        event_file = ['']
        
    return event_file[0]

def covariate_frame(data):
    """Grab dataframe containing covariate information"""
    from bids.layout import BIDSLayout
    layout = BIDSLayout(data)
    file = layout.get(return_type='filename', extension='.tsv', suffix='participants')
    
    return file[0]

def split_(smoothed, unsmoothed, half):
    """Split functional scan into first and second half for NPAIRs genetic algorithm"""
    import nibabel as nib
    from nipype.interfaces.fsl import ExtractROI
    length = round(nib.load(smoothed).shape[-1] / 2)
    
    if half == 'first':
        smoothed = ExtractROI(in_file=smoothed, t_min=0, t_size=length).run().outputs.roi_file
        unsmoothed = ExtractROI(in_file=unsmoothed, t_min=0, t_size=length).run().outputs.roi_file
    elif half == 'second':
        smoothed = ExtractROI(in_file=smoothed, t_min=length, t_size=-1).run().outputs.roi_file
        unsmoothed = ExtractROI(in_file=unsmoothed, t_min=length, t_size=-1).run().outputs.roi_file
    else:
        raise ValueError("Only 'first' and 'second' are valid inputs for split half analysis, but {half} was given instead.".format(half=half))
        
    return smoothed, unsmoothed
    
        
def remove(T1w, bold):
    """Remove container"""
    return T1w[0], bold[0]
    
def insert(string, ind, new):
    """Insert string at specified location"""
    return string[:ind] + new + string[ind:]

def make_buff_vars(dic):
    """Creates dynamic functions to index parameters for the correct pipelines - allows for shared data as long as possible"""
    func = "def buff_var({var}):\n\treturn "
    var = [param_key for key in dic for param_key in dic[key] if param_key != 'id' or param_key != 'link']
    inputs = ''
    ret = ''
    for v in var:
        inputs += str(v) + ', '
        ret += str(v) + '[i-i_sub], '
        
    inputs += 'i, i_in, i_sub'
    ret += 'i'
    input_names = inputs.split(', ')
    return func.format(var=inputs) + ret, input_names, input_names[:-2]
        
def setatts(node, dic, keys):
    """Sets the inputs for automated buffer nodes"""
    keys = keys[:-3]
    for key in keys:
        setattr(node.inputs, key, dic[key])
        
def get_links(dic, keys):
    """Creates dictionary of links between pipelines (i.e. when a pipeline splits off from its progenitor)"""
    nodes = []
    connections = {}
    expansion = {}
    for key in keys:
        check = False
        sub_dic = dic[key]
        try:
            #check if pipeline is the parent of any future pipelines
            if len(sub_dic['id']) > 1:
                check = True
                pipeline_keys = list(sub_dic.keys())
                pipeline_keys.remove('id')
            else:
                pipeline_keys = list(sub_dic[key].keys())
                pipeline_keys.remove('id')
                #no dependent pipelines, add link directly
                if pipeline_keys[0] not in connections:
                    connections[pipeline_keys[0]] = {}
                connections[pipeline_keys[0]][key] = [key]
        except:
            #id missing, or dic[key][key] doesn't exist
            check = True
            pipeline_keys = sorted(list(sub_dic.keys()))
            if pipeline_keys == ['id']:
                pipeline_keys = []
        #make connection links
        if check:
            for k in pipeline_keys:
                #splits off from parent pipeline
                if 'link' in sub_dic[k]:
                    link = sub_dic[k]['link']
                    node = list(link.values())[0][0]
                    link_key = list(link.keys())[0]
                    key_ = link_key.copy()
                    nodes.append(node)
                    try:
                        if link_key in connections[node]:
                            connections[node][link_key].append(k)
                        else:
                            current_vals = connections[node].values()
                            current_keys = list(connections[node].keys())
                            #check to ensure link to root pipeline, not intermediate (i.e. parameters different at later stage, but node the same)
                            for i, vals in enumerate(current_vals):
                                if link_key in vals:
                                    key_ = current_keys[i]
                                    break
                            
                            if key_ in connections[node]:
                                connections[node][key_].append(k)
                            else:
                                connections[node][key_] = [link_key]
                                #ensure pipelines link to themselves
                                if link_key != k:
                                    connections[node][key_].append(k)
                    except:
                        #either node not in connections, or link_key not in connections[node]
                        if node in connections:
                            connections[node][link_key] =  [link_key]
                        else:
                            connections[node] = {link_key: [link_key]}
                        
                        if k not in connections[node][link_key]:
                            connections[node][link_key].append(k)
                else:
                    #link pipeline to itself
                    pipeline_keys = list(sub_dic[key].keys())
                    pipeline_keys.remove('id')
                    if pipeline_keys[0] not in connections:
                        connections[pipeline_keys[0]] = {}
                    connections[pipeline_keys[0]][key] = [key]
    
    return list(dict.fromkeys(nodes)), add_mapping(connections)

def add_mapping(con):
    """Ensure that pipelines with no children are included (self link)"""
    single = []
    out_con = con.copy()
    for key in con:
        values = [item for value in con[key].values() for item in value]
        if single:
            missing = sorted((Counter(single) - Counter(values)).elements())
            for pipe in missing:
                out_con[key][pipe] = [pipe]
            single = [item for value in out_con[key].values() for item in value]
        else:
            single = values
    return out_con


def traverse(dic, flow, suffix, pipeline, to_run):
    """Makes connections between nodes, sets changeable parameters of nodes"""
    #check to see if there are repeat pipelines to be avoided, get starting pipeline
    if to_run:
        start_pipe_ = min(to_run)
    else:
        start_pipe_ = pipeline
    #Node to start process of differentiation into unique pipelines (everything splits off from start_pipe_ pipeline)
    iternode = Node(IdentityInterface(fields=['i']), name='iternode'+suffix)
    iternode.iterables = ('i', [start_pipe_])
    
    buff_count = []
    #iterate through workflow stages (i.e. preprocessing, level 1, etc.)
    for wf in dic:
        buff_dic = {}
        dic_ = dic[wf]
        dic_k = list(dic_.keys())
        start_pipe = [i for i in dic_k if isinstance(i, (int, np.integer))]
        split_nodes, connections = get_links(dic_, start_pipe)
        outstanding = False
        #check if pipelines split in wf stage
        if connections:
            for i, info in enumerate(dic_[start_pipe_][start_pipe_]):
                if info == 'id' or info == 'link':
                    continue
                #check if first node, or if node is where pipelines split
                if info not in split_nodes or not i:
                    buff_dic[info] = dic_[start_pipe_][start_pipe_][info]
                    outstanding = True
                else:
                    #create a buffer function to index parameters for pipelines
                    if not buff_count:
                        buff_count.append(1)
                    connections[buff_count[-1]] = connections.pop(list(buff_dic.keys())[0])
                    func, input_names, output_names = make_buff_vars(buff_dic)
                    vars()['buff_' + str(buff_count[-1])] = Node(Function(input_names=input_names, output_names=output_names), name='buff_' + str(buff_count[-1]))
                    vars()['buff_' + str(buff_count[-1])].inputs.function_str = func
                    vars()['buff_' + str(buff_count[-1])].inputs.i_sub = start_pipe_
                    setatts(vars()['buff_' + str(buff_count[-1])], dic_, input_names)
                    #make connections between buffer node and parameters of destination nodes
                    for name in input_names[:-3]:
                        #nodename_parametername
                        end = re.search('^([A-Za-z0-9]+)_([A-Za-z_]+)', name)
                        if flow.get_node(wf).get_node(end.group(1)):
                            flow.get_node(wf).connect(vars()['buff_' + str(buff_count[-1])], name, flow.get_node(wf).get_node(end.group(1)), end.group(2))
                    buff_count.append(buff_count[-1] + 1)
                    #reset buff_dic with new starting node information (split node)
                    buff_dic = {info: dic_[start_pipe_][start_pipe_][info]}
                    outstanding = True
        #check if there are connections which haven't been made
        if outstanding:
            if not buff_count:
                buff_count.append(1)
            #create buffer function
            connections[buff_count[-1]] = connections.pop(list(buff_dic.keys())[0])
            func, input_names, output_names = make_buff_vars(buff_dic)
            vars()['buff_' + str(buff_count[-1])] = Node(Function(input_names=input_names, output_names=output_names), name='buff_' + str(buff_count[-1]))
            vars()['buff_' + str(buff_count[-1])].inputs.function_str = func
            vars()['buff_' + str(buff_count[-1])].inputs.i_sub = start_pipe_
            setatts(vars()['buff_' + str(buff_count[-1])], dic_, input_names)
            #make connections
            for name in input_names[:-3]:
                end = re.search('^([A-Za-z0-9]+)_([A-Za-z_]+)', name)
                if flow.get_node(wf).get_node(end.group(1)):
                    flow.get_node(wf).connect(vars()['buff_' + str(buff_count[-1])], name, flow.get_node(wf).get_node(end.group(1)), end.group(2))
            #connect i to all buff functions -> used to index parameters as defined with connections
            for buff in buff_count:
                if buff == 1 and vars().get('buff_1', False):
                    vars()['buff_' + str(buff)].itersource = ('iternode'+suffix, 'i')
                    vars()['buff_' + str(buff)].iterables = [('i', connections[buff])]
                    flow.get_node(wf).connect(iternode, 'i', vars()['buff_' + str(buff)], 'i_in')
                elif not vars().get('buff_1', False):
                    break
                else:
                    vars()['buff_' + str(buff)].itersource = ('buff_' + str(buff - 1), 'i')
                    vars()['buff_' + str(buff)].iterables = [('i', connections[buff])]
                    flow.get_node(wf).connect(vars()['buff_' + str(buff - 1)], 'i', vars()['buff_' + str(buff)], 'i_in')
        #Connect all parameters that are configurable, but don't change
        for node in dic_['const']:
            const = dic_['const'][node]
            keys = list(const.keys())
            vals = list(const.values())
            for i, k in enumerate(keys):
                k_var = re.search(node+'_([A-Za-z_]+)', k).group(1)
                
                if flow.get_node(wf).get_node(node) == None:
                    continue
                
                if k_var in flow.get_node(wf).get_node(node).inputs.get():
                    setattr(flow.get_node(wf).get_node(node).inputs, k_var, vals[i])
        
        if buff_count and outstanding:
            buff_count = [buff_count[-1]+1]
                

def define_paths(container, dictionary, indexes):
    """Formats dictionary used to store parameters and keeps track of where pipelines split"""
    out_dic = {}
    link = {}
    
    #commented block should be ok to remove (doesn't seem to do anything), but haven't tested so left for now
    #old_x = len(indexes)
    #get starting index of last pipeline
    #if old_x >= 1:
    #    old_x -= 1
        
    #change = len(indexes) - 1
    
    #initialize dictionary
    for i, vals in enumerate(indexes):
        if isinstance(vals, np.ndarray):
            out_dic[int(indexes[i].min())] = {'id': vals}
            rng = indexes[0]
        else:
            if type(vals) == tuple:
                rng = np.array(range(vals[0], vals[1]))
            else:
                rng = np.array(range(vals))
                
            out_dic[rng[i]] = {'id': rng}
            
    out_dic['const'] = {}
    #start pipeline
    for i, key in enumerate(dictionary):
        #parameters
        for subkey in dictionary[key]:
            try:
                #change every element to string so that values can be sorted
                placeholder = [str(element) for element in dictionary[key][subkey]]
                container = np.vstack((container, placeholder))
            except:
                pass
            #Find pipelines which are identical
            vals, ind = np.unique(container, return_inverse=True, axis=1)
            index = [np.where((vals[:,i].reshape(-1,1) == container).sum(axis=0) == container.shape[0])[0] for i in range(vals.shape[1])]
            index_ = sorted(index, key=min)
            start_ = rng[0]
            #adjust indices to be global, not just for iteration of GA
            index_ = [arr + start_ for arr in index_]
            #convert datatype if necessary to check if constant across pipelines
            try:
                gen = np.unique(dictionary[key][subkey])
            except:
                tostring = [str(item) for item in dictionary[key][subkey]]
                gen = np.unique(tostring)
            #iterate through root/starting pipelines for each stage
            for k in out_dic:
                if not isinstance(k, (int, np.integer)):
                    break
                #add constant
                if len(gen) == 1:
                    if key not in out_dic['const']:
                        out_dic['const'][key] = {}
                    if subkey not in out_dic['const'][key]:
                        out_dic['const'][key][subkey] = {}
                                        
                    out_dic['const'][key][subkey] = dictionary[key][subkey][0]
                    continue
                #iterate through pipelines
                for j, x in enumerate(index_):
                    #don't classify pipelines if they have yet to split off
                    if 'id' in out_dic[k]:
                        if len(np.intersect1d(x, out_dic[k]['id'])) != len(x):
                            continue
                    #safely construct dictionary
                    if min(index_[j]) not in out_dic[k]:
                        out_dic[k][min(index_[j])] = {}
                    if key not in out_dic[k][min(index_[j])]:
                        out_dic[k][min(index_[j])][key] = {}
                    if subkey not in out_dic[k][min(index_[j])][key]:
                        out_dic[k][min(index_[j])][key][subkey] = {}
                    #flags as initial run
                    if 'id' not in out_dic[k][min(index_[j])]:
                        out_dic[k][min(index_[j])]['id'] = [-1]
                    #special case for if pipeline has not been added
                    #create directory of links between pipelines
                    if not np.array_equiv(out_dic[k][min(index_[j])]['id'], x):
                        #inital run check
                        if not np.array_equiv(out_dic[k][min(index_[j])]['id'], [-1]):
                            cx = Counter(out_dic[k][min(index_[j])]['id'])
                            cid = Counter(x)
                            #create links to parent pipeline
                            for out in sorted((cx - cid).elements()):
                                link[out] = {min(index_[j]): [key, subkey]}
                        #create link if splits off at intermediate stage
                        elif 'id' in out_dic[k]:
                            if not np.array_equiv(out_dic[k][min(index_[j])]['id'], [-1]) or (min(x) == min(out_dic[k]['id']) and len(x) < len(out_dic[k]['id'])):
                                for out in out_dic[k]['id']:
                                    link[out] = {min(index_[j]): [key, subkey]}
                    #save link in dictionary
                    if min(index_[j]) in link:
                        out_dic[k][min(index_[j])]['link'] = link[min(index_[j])]
                        
                    out_dic[k][min(index_[j])]['id'] = x
                    out_dic[k][min(index_[j])][key][subkey] = None
                #Below commented out, am pretty sure can be deleted, but left in case I'm wrong
                #if k == change:
                #    old_x = j
                                
            link = {}
            out_dic[subkey] = dictionary[key][subkey]
        
    return container, out_dic, index_

def load(path, file):
    out = os.path.join(out_dir, path, file)
    if os.path.isfile(out):
        with open(out, 'rb') as f:
            loaded = pickle.load(f)
    else:
        loaded = ''
    
    return loaded

def save(path, file, frame):
    out = os.path.join(out_dir, path)
    Path(out).mkdir(parents=True, exist_ok=True)
    out = os.path.join(out, file)
        
    with open(out, 'wb') as f:
        pickle.dump(frame, f)
    
    return out

def organize(task, out_frame):
    """Creates a dictionary of final output files, and parameters for each pipeline - excludes parameters that are unchanged across all pipelines
       
       Structure:
           {pipeline: {network: {contrast: file}},
                      {parameters: {parameters}}
                      }
    """
    processed = {'pipeline': {}}
    pathlist = Path(out_dir+'/pipelines/'+task).glob('**/*_corrected_[0-9]*')
    dat_frame = out_frame
    
    with open(dat_frame, 'rb') as file:
        dat_frame = pickle.load(file)
        
    comp = pd.DataFrame(np.roll(dat_frame.values, 1, axis=0), index=dat_frame.index)
        
    for path in pathlist:
        path = str(path)
        network = int(re.search('.*_network_([0-9]+)', path).group(1))
        contrast = int(re.search('.*_corrected_([0-9]+).nii.gz', path).group(1))
        pipeline = int(re.search('.*_i_([0-9]+)', path).group(1))
        if pipeline in processed['pipeline']:
            if network in processed['pipeline'][pipeline]['network']:
                processed['pipeline'][pipeline]['network'][network]['contrast'][contrast] = path
            else:
                processed['pipeline'][pipeline]['network'][network] = {'contrast': {contrast: path}}
        else:
            processed['pipeline'][pipeline] = {'network': {network: {'contrast': {contrast: path}}}}
            
        pipe_dat = dat_frame.loc[pipeline]
        for i, column in enumerate(dat_frame):
            col = pipe_dat[column]
            if (comp[i] == dat_frame[column]).all():
                continue
            
            if 'parameters' not in processed['pipeline'][pipeline]:
                processed['pipeline'][pipeline]['parameters'] = {}
            
            if isinstance(col, dict):
                for key in col:
                    processed['pipeline'][pipeline]['parameters'][key] = col[key]
            else:
                processed['pipeline'][pipeline]['parameters'][column] = col
    
    return save('', task+'_organized.pkl', processed)











