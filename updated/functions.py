#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:54:23 2021

@author: grahamseasons
"""
from nipype import Node, IdentityInterface, Function
import re
from collections import Counter
import numpy as np

def make_buff_vars(dic):
    func = "def buff_var({var}):\n\treturn "
    var = [param_key for key in dic for param_key in dic[key] if param_key != 'id' or param_key != 'link']
    inputs = ''
    ret = ''
    for v in var:
        inputs += str(v) + ', '
        ret += str(v) + '[i], '
        
    inputs += 'i, i_in' #inputs.rstrip(', ')
    ret += 'i' #ret.rstrip(', ')
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
                #pipeline_keys.sort()
        except:
            check = True
            pipeline_keys = sorted(list(sub_dic.keys()))
        
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
                                connections[node][key_].append(k)
                            #index = list(connections[list(link.values())[0][0]].keys())[0]
                            #connections[list(link.values())[0][0]][index].append(k)
                    except:
                        if node in connections:
                            connections[node][link_key] =  [link_key]
                        else:
                            connections[node] = {link_key: [link_key]}
                        
# =============================================================================
#                         if list(link.keys())[0] in list(connections[list(link.values())[0][0]].values()):
#                             index = list(connections[list(link.values())[0][0]].keys())[0]
#                             connections[list(link.values())[0][0]][index].append(k)
# =============================================================================
                        
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
        
#will have suffix definition for iternode [x]
#change so that every connection finds node within flow inputs [x]
#change so that constant nodes are added to workflow before being passed to traverse []
#still need to add buff for inputs to workflow(up a level) and connect iternode to it []
#no underscores in node names
def traverse(dic, flow, suffix):
    dic_k = list(dic.keys())
    start_pipe = [i for i in dic_k if type(i) == int]
    split_nodes, connections = get_links(dic, start_pipe)
    iternode = Node(IdentityInterface(fields=['i']), name='iternode'+suffix)
    iternode.iterables = ('i', start_pipe)
    
    buff_count = []
    buff_dic = {}
    for i, info in enumerate(dic[0][0]):
        if info == 'id' or info == 'link':
            continue
        if info not in split_nodes or not i:
            buff_dic[info] = dic[0][0][info]
            outstanding = True
        else:
            if not buff_count:
                buff_count.append(1)
            connections[buff_count[-1]] = connections.pop(list(buff_dic.keys())[0])
            func, input_names, output_names = make_buff_vars(buff_dic)
            vars()['buff_' + str(buff_count[-1])] = Node(Function(input_names=input_names, output_names=output_names), name='buff_' + str(buff_count[-1]))
            vars()['buff_' + str(buff_count[-1])].inputs.function_str = func
            setatts(vars()['buff_' + str(buff_count[-1])], dic, input_names)
            for name in input_names[:-2]:
                end = re.search('^([A-Za-z]+)_([A-Za-z_]+)', name)
                flow.connect(vars()['buff_' + str(buff_count[-1])], name, flow.get_node(end.group(1)), end.group(2))
            buff_count.append(buff_count[-1] + 1)
            buff_dic = {info: dic[0][0][info]}
            outstanding = False
            
    if outstanding:
        if not buff_count:
            buff_count.append(1)
        connections[buff_count[-1]] = connections.pop(list(buff_dic.keys())[0])
        func, input_names, output_names = make_buff_vars(buff_dic)
        vars()['buff_' + str(buff_count[-1])] = Node(Function(input_names=input_names, output_names=output_names), name='buff_' + str(buff_count[-1]))
        vars()['buff_' + str(buff_count[-1])].inputs.function_str = func
        setatts(vars()['buff_' + str(buff_count[-1])], dic, input_names)
        #buff_count.append(buff_count[-1] + 1)
    
        for name in input_names[:-2]:
            end = re.search('^([A-Za-z]+)_([A-Za-z_]+)', name)
            flow.connect(vars()['buff_' + str(buff_count[-1])], name, flow.get_node(end.group(1)), end.group(2))
    
    for buff in buff_count:
        if buff == 1:
            vars()['buff_' + str(buff)].itersource = ('iternode'+suffix, 'i')
            vars()['buff_' + str(buff)].iterables = [('i', connections[buff])]
            flow.connect(iternode, 'i', vars()['buff_' + str(buff)], 'i_in')
        else:
            vars()['buff_' + str(buff)].itersource = ('buff_' + str(buff - 1), 'i')
            vars()['buff_' + str(buff)].iterables = [('i', connections[buff])]
            flow.connect(vars()['buff_' + str(buff - 1)], 'i', vars()['buff_' + str(buff)], 'i_in')
    
    for node in dic['const']:
        const = dic['const'][node]
        keys = list(const.keys())
        vals = list(const.values())
        for i, k in enumerate(keys):
            k_var = re.search(node+'_([A-Za-z_]+)', k).group(1)
            setattr(flow.get_node(node).inputs, k_var, vals[i])
    
    if buff_count:
        joinsource = 'buff_' + str(buff)
    else:
        joinsource = ''
        
    return joinsource


def format_out(smoothed, segmentations, warp_file, outliers, mc_par, brain):
    eff=3


def define_paths(container, dictionary, indexes):
    out_dic = {}
    link = {}
    old_x = len(indexes)
    
    if old_x >= 1:
        old_x -= 1
        
    change = len(indexes) - 1
    
    for i, vals in enumerate(indexes):
        if isinstance(vals, np.ndarray):
            out_dic[i] = {'id': vals}
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
                
            gen = np.unique(dictionary[key][subkey])
            
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
                                
# =============================================================================
#                     if j > old_x:
#                         try:
#                             out_dic[k][min(index_[j])]['link'] = link[min(index_[j])]
#                         except:
#                             A=3
# =============================================================================
                        
                    out_dic[k][min(index_[j])]['id'] = x
                    out_dic[k][min(index_[j])][key][subkey] = None
                        
                if k == change:
                    old_x = j
                                
            link = {}
            out_dic[subkey] = dictionary[key][subkey]
        
    return container, out_dic, index_















