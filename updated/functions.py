#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:54:23 2021

@author: grahamseasons
"""
from nipype import Node, IdentityInterface, Function
from nipype.interfaces.base import Undefined
import re, random, os
from collections import Counter
import numpy as np

def no_mask(file):
    raise FileNotFoundError("Specified mask '{m}' does not exist".format(m=file))
    
def invalid(mapped):
    raise SyntaxError("Input paramater '{name}' in an unsupported format. Acceptable formats are (brain region, thr), (brain region, thr_min, thr_max, or path_to_mask".format(name=mapped))

def generate_dictionaries(map_genes, links, params, pop, multiscan, wiggle):
    container = np.zeros((1, pop.shape[0]), str)
    indexes = [pop.shape[0]]
    expand_inputs = {}
    
    preprocess = {}
    level1 = {}
    
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
            for link in links[dic]:
                if previous in link[:len(previous)]:
                    connect = links[dic][link]
                    if type(connect) == list:
                        group = re.search('([A-Za-z0-9]+)_', connect[0]).group(1)
                        if group not in vars()[dic]:
                            for opt in valid:
                                if group in vars()[opt]:
                                    break
                        vals = vars()[dic][group][connect[0]]
                        check = vars()[dic][group][connect[1]]
                        if len(connect) == 3:
                            rule = connect[2]
                            vars()[dic][previous][link] = [val if check[c] != rule else Undefined for c, val in enumerate(vals)]
                        else:
                            vars()[dic][previous][link] = [val if check[c] else Undefined for c, val in enumerate(vals)]
                    else:
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
                            
            if 'F' == previous[0]:
                expand_inputs[previous[1:]] = vars()[dic][previous]
                
            previous = node_name
            
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
        
        if node_name not in vars()[dic]:
            vars()[dic][node_name] = {}
            
        values = params[key-counter,:]
        
        isint = False
        if values.dtype == float:
            check_vals = [val.is_integer() for val in values]
            if sum(check_vals) == len(check_vals):
                isint = True
        
        if keys[0][0] == '!':
            continue
        
        if len(gene) > 1 or '~construct~' in keys[0]:
            for l, i in enumerate(values):
                if round(i) in gene:
                    mapped = gene[round(i)]
                else:
                    mapped = i
                if type(mapped) == dict and keys[0][-1] == '_':
                    for k in mapped:
                        if keys[0][-1] == '_':
                            param = keys[0] + k
                        else:
                            param = node_name + '_' + k
                        if param not in vars()[dic][node_name]:
                             vars()[dic][node_name][param] = []
                             
                        vars()[dic][node_name][param].append(mapped[k])
                elif '~construct~' in keys[0]:
                    var_name = re.search('_([A-Za-z]+)', keys[0]).group(1)
                    key_name = re.search('_([A-Za-z]+)$', keys[0]).group(1)
                    param = node_name + '_' + var_name
                    
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    
                    if not isint:
                        rand = random.Random(i)
                        if type(mapped) == list:
                            mapped = [(m[0], rand.randint(m[1]-wiggle if m[1]-wiggle > 0 else 0, m[1]+wiggle if m[1]+wiggle > 95 else 95)) if len(m) == 2 and isinstance(m, (tuple,list)) else (m[0], rand.randint(m[1], m[2])) if isinstance(m, (tuple,list)) else m if isinstance(m, str) else invalid(param+'_'+key_name) for m in mapped]
# =============================================================================
#                             if type(mapped[0]) == tuple or type(mapped[0]) == list:
#                                 mapped = [(m[0], rand.randint(m[1]-wiggle if m[1]-wiggle > 0 else 0, m[1]+wiggle if m[1]+wiggle > 95 else 95)) if len(m) == 2 else (m[0], rand.randint(m[1], m[2])) for m in mapped]
#                             else:
#                                 mapped = [file if os.path.isfile(file) else no_mask(file) for file in mapped]
# =============================================================================
                        elif mapped == tuple:
                            mapped = rand.randint(mapped[0], mapped[1])
                        else:
                            mapped = i
                    else:
                        if isinstance(mapped, float):
                            mapped = int(mapped)
                            
                    if len(vars()[dic][node_name][param]) != pop.shape[0]:
                        construction_key = key_name
                        vars()[dic][node_name][param].append({key_name: mapped})
                    else:
                        if key_name in links[dic]:
                            if vars()[dic][node_name][param][l][construction_key] not in links[dic][key_name]:
                                continue
                            
                        vars()[dic][node_name][param][l][key_name] = mapped
                else:
                    param = keys[0]
                        
                    if param not in vars()[dic][node_name]:
                        vars()[dic][node_name][param] = []
                    
                    vars()[dic][node_name][param].append(mapped)
        else:
            if values.dtype == float:
                if isint:
                    values = values.astype(int)
                    
            vars()[dic][node_name][keys[0]] = values
            
    return master, expand_inputs
    
def metadata(filename, data):
    from bids.layout import BIDSLayout
    layout = BIDSLayout(data)
    meta_data = layout.get_metadata(filename)
    TR = meta_data['RepetitionTime']
    return TR

def event_grabber(file, data):
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
    from bids.layout import BIDSLayout
    layout = BIDSLayout(data)
    file = layout.get(return_type='filename', extension='.tsv', suffix='participants')
    
    return file[0]

def split_(smoothed, unsmoothed, half):
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
    return T1w[0], bold[0]
    
def insert(string, ind, new):
    return string[:ind] + new + string[ind:]

def make_buff_vars(dic):
    func = "def buff_var({var}):\n\treturn "
    var = [param_key for key in dic for param_key in dic[key] if param_key != 'id' or param_key != 'link']
    inputs = ''
    ret = ''
    for v in var:
        inputs += str(v) + ', '
        ret += str(v) + '[i], '
        
    inputs += 'i, i_in'
    ret += 'i'
    input_names = inputs.split(', ')
    return func.format(var=inputs) + ret, input_names, input_names[:-1]
        
def setatts(node, dic, keys):
    keys = keys[:-2]
    for key in keys:
        setattr(node.inputs, key, dic[key])
        
def get_links(dic, keys):
    nodes = []
    connections = {}
    for key in keys:
        check = False
        sub_dic = dic[key]
        try:
            if len(sub_dic['id']) > 1:
                check = True
                pipeline_keys = list(sub_dic.keys())
                pipeline_keys.remove('id')
            else:
                pipeline_keys = list(sub_dic[key].keys())
                if pipeline_keys[0] not in connections:
                    connections[pipeline_keys[0]] = {}
                connections[pipeline_keys[0]][key] = [key]
        except:
            check = True
            pipeline_keys = sorted(list(sub_dic.keys()))
            if pipeline_keys == ['id']:
                pipeline_keys = []
        
        if check:
            for k in pipeline_keys:
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
                            for i, vals in enumerate(current_vals):
                                if link_key in vals:
                                    key_ = current_keys[i]
                                    break
                            
                            if key_ in connections[node]:
                                connections[node][key_].append(k)
                            else:
                                connections[node][key_] = [link_key]
                                if link_key != k:
                                    connections[node][key_].append(k)
                    except:
                        if node in connections:
                            connections[node][link_key] =  [link_key]
                        else:
                            connections[node] = {link_key: [link_key]}
                        
                        if k not in connections[node][link_key]:
                            connections[node][link_key].append(k)
    
    return list(dict.fromkeys(nodes)), add_mapping(connections)

def add_mapping(con):
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


def traverse(dic, flow, suffix):
    start_pipe = 0
    iternode = Node(IdentityInterface(fields=['i']), name='iternode'+suffix)
    iternode.iterables = ('i', [start_pipe])
    
    buff_count = []
    for wf in dic:
        buff_dic = {}
        dic_ = dic[wf]
        dic_k = list(dic_.keys())
        start_pipe = [i for i in dic_k if type(i) == int]
        split_nodes, connections = get_links(dic_, start_pipe)
        outstanding = False
        if connections:
            for i, info in enumerate(dic_[0][0]):
                if info == 'id' or info == 'link':
                    continue
                if info not in split_nodes or not i:
                    buff_dic[info] = dic_[0][0][info]
                    outstanding = True
                else:
                    if not buff_count:
                        buff_count.append(1)
                    connections[buff_count[-1]] = connections.pop(list(buff_dic.keys())[0])
                    func, input_names, output_names = make_buff_vars(buff_dic)
                    vars()['buff_' + str(buff_count[-1])] = Node(Function(input_names=input_names, output_names=output_names), name='buff_' + str(buff_count[-1]))
                    vars()['buff_' + str(buff_count[-1])].inputs.function_str = func
                    setatts(vars()['buff_' + str(buff_count[-1])], dic_, input_names)
                    for name in input_names[:-2]:
                        end = re.search('^([A-Za-z0-9]+)_([A-Za-z_]+)', name)
                        flow.get_node(wf).connect(vars()['buff_' + str(buff_count[-1])], name, flow.get_node(wf).get_node(end.group(1)), end.group(2))
                    buff_count.append(buff_count[-1] + 1)
                    buff_dic = {info: dic_[0][0][info]}
                    outstanding = True
                
        if outstanding:
            if not buff_count:
                buff_count.append(1)
            connections[buff_count[-1]] = connections.pop(list(buff_dic.keys())[0])
            func, input_names, output_names = make_buff_vars(buff_dic)
            vars()['buff_' + str(buff_count[-1])] = Node(Function(input_names=input_names, output_names=output_names), name='buff_' + str(buff_count[-1]))
            vars()['buff_' + str(buff_count[-1])].inputs.function_str = func
            setatts(vars()['buff_' + str(buff_count[-1])], dic_, input_names)
        
            for name in input_names[:-2]:
                end = re.search('^([A-Za-z0-9]+)_([A-Za-z_]+)', name)
                flow.get_node(wf).connect(vars()['buff_' + str(buff_count[-1])], name, flow.get_node(wf).get_node(end.group(1)), end.group(2))
        
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
        
        for node in dic_['const']:
            const = dic_['const'][node]
            keys = list(const.keys())
            vals = list(const.values())
            for i, k in enumerate(keys):
                k_var = re.search(node+'_([A-Za-z_]+)', k).group(1)
                setattr(flow.get_node(wf).get_node(node).inputs, k_var, vals[i])
        
        if buff_count and outstanding:
            buff_count = [buff_count[-1]+1]
                

def define_paths(container, dictionary, indexes):
    out_dic = {}
    link = {}
    old_x = len(indexes)
    
    if old_x >= 1:
        old_x -= 1
        
    change = len(indexes) - 1
    
    for i, vals in enumerate(indexes):
        if isinstance(vals, np.ndarray):
            out_dic[int(indexes[i].min())] = {'id': vals}
        else:
            out_dic[i] = {'id': np.array(range(vals))}
            
    out_dic['const'] = {}
    
    for i, key in enumerate(dictionary):
        for subkey in dictionary[key]:
            try:
                placeholder = [str(element) for element in dictionary[key][subkey]]
                container = np.vstack((container, placeholder))
            except:
                A=3
            vals, ind = np.unique(container, return_inverse=True, axis=1)
            index = [np.where((vals[:,i].reshape(-1,1) == container).sum(axis=0) == container.shape[0])[0] for i in range(vals.shape[1])]
            index_ = sorted(index, key=min)
               
            try:
                gen = np.unique(dictionary[key][subkey])
            except:
                tostring = [str(item) for item in dictionary[key][subkey]]
                gen = np.unique(tostring)
            
            for k in out_dic:
                if type(k) != int:
                    break
                if len(gen) == 1:
                    if key not in out_dic['const']:
                        out_dic['const'][key] = {}
                    if subkey not in out_dic['const'][key]:
                        out_dic['const'][key][subkey] = {}
                                        
                    out_dic['const'][key][subkey] = dictionary[key][subkey][0]
                    continue
                                    
                for j, x in enumerate(index_):
                    if 'id' in out_dic[k]:
                        if len(np.intersect1d(x, out_dic[k]['id'])) != len(x):
                            continue
                    if min(index_[j]) not in out_dic[k]:
                        out_dic[k][min(index_[j])] = {}
                    if key not in out_dic[k][min(index_[j])]:
                        out_dic[k][min(index_[j])][key] = {}
                    if subkey not in out_dic[k][min(index_[j])][key]:
                        out_dic[k][min(index_[j])][key][subkey] = {}
                    if 'id' not in out_dic[k][min(index_[j])]:
                        out_dic[k][min(index_[j])]['id'] = [-1]
                        
                    if not np.array_equiv(out_dic[k][min(index_[j])]['id'], x):
                        if not np.array_equiv(out_dic[k][min(index_[j])]['id'], [-1]):
                            cx = Counter(out_dic[k][min(index_[j])]['id'])
                            cid = Counter(x)
                            if len(sorted((cx - cid).elements())) > 1:
                                A=3
                            for out in sorted((cx - cid).elements()):
                                link[out] = {min(index_[j]): [key, subkey]}
                        elif 'id' in out_dic[k]:
                            if not np.array_equiv(out_dic[k][min(index_[j])]['id'], [-1]) or (min(x) == min(out_dic[k]['id']) and len(x) < len(out_dic[k]['id'])):
                                for out in out_dic[k]['id']:
                                    link[out] = {min(index_[j]): [key, subkey]}
                                    
                    if min(index_[j]) in link:
                        out_dic[k][min(index_[j])]['link'] = link[min(index_[j])]
                        
                    out_dic[k][min(index_[j])]['id'] = x
                    out_dic[k][min(index_[j])][key][subkey] = None
                        
                if k == change:
                    old_x = j
                                
            link = {}
            out_dic[subkey] = dictionary[key][subkey]
        
    return container, out_dic, index_















