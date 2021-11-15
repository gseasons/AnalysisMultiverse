#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:45:59 2021

@author: grahamseasons
"""
# =============================================================================
import subprocess
import os, sys, re, glob
from shutil import copy2
import xml.etree.ElementTree as ET
import argparse
import json
#from updated.parameters import genes
# from kivy.app import App
# from kivy.uix.label import Label
# 
# class multiverse(App):
#     def build(self):
#         root_widget = Label(text='')
#         return root_widget
#     
# multiverse().run()
# =============================================================================
data_dir = '/Volumes/NewVolume/super_agers'
exp_dir = '/Volumes/NewVolume/sup_pre_full'
working_dir = 'working_dir'
out_dir = exp_dir + '/processed'
code = '/Users/grahamseasons/fMRI/analysis_multiverse/updated'

import tkinter as tk
from tkinter import ttk
from functools import partial

class MultiverseConfig:
    def __init__(self):
        self.master = tk.Toplevel()
        self.out_dic = {'level1': {'~construct~Finfo_rest_seedinfo': {'genes': []}, '~construct~Finfo_rest_type': {'genes': [], 0: 'atlas', 1: 'data', 2: 'ROI'}}}
        self.config_dic = {}
        self.tabcontrol = ttk.Notebook(self.master)
        self.tab_main = ttk.Frame(self.tabcontrol)
        self.tab_main.columnconfigure(0, weight=1)
        self.tabcontrol.add(self.tab_main, text='Essential')
        
        self.essential()
        label_dic = {'preprocess': 'Preprocessing', 'level1': 'Level 1 Analysis'}
        data = self.get_format()
        for i, stage in enumerate(data):
            counter = 0
            self.out_dic[stage] = {}
            vars(self)['tab_'+stage] = ttk.Frame(self.tabcontrol)
            vars(self)['tab_'+stage].columnconfigure(0, weight=1)
            self.tabcontrol.add(vars(self)['tab_'+stage], text=label_dic[stage])
            vars(self)['tab_window_'+stage] = ttk.Frame(vars(self)['tab_'+stage])
            vars(self)['tab_window_'+stage].grid(row=0, column=0, sticky="nsew")
            vars(self)['tab_window_'+stage].columnconfigure(0, weight=1)
            counter += 1
            for node in data[stage]:
                ttk.Label(vars(self)['tab_window_'+stage], text=node['alias'], font='Helvetica 12 bold').grid(row=counter)
                counter += 1
                for param in node['params']:
                    param_name = param['name']
                    alias = param['alias']
                    default = param['default']
                    node_name = node['node_name']
                    param_name_full = node_name + '_' + param_name
                    on_off = param['on_off']
                    
                    if 'value' in param:
                        value = param['value']
                        value_map = ''
                    else:
                        value_map = param['value_map']
                        value = list(range(len(value_map)))
                        
                    if on_off:
                        butt_text = 'On'
                        self.out_dic[stage][param_name_full] = {'gene': value}
                    else:
                        butt_text = 'Off'
                        self.out_dic[stage][param_name_full] = {'gene': value[default]}
                        
                    if value_map:
                        for i in value:
                            val = value_map[i]
                            if type(val) == list:
                                val = tuple(val)
                            self.out_dic[stage][param_name_full][i] = val
                            
                        self.strings(stage, alias, default, value_map, param_name_full, value, butt_text, counter)
                    else:
                        self.numerical(stage, alias, default, value_map, param_name_full, value, butt_text, counter)
                    
                    counter += 6
                break
            break

        self.tabcontrol.pack(expand=1, fill="both")
        tk.mainloop()
        A = 3 #SET UP DEFUALT VALUES IN DICTIONARY AS CREATING CONFIG BUTTONS
        #ADD SEEDS ADDED OPTION FOR EACH METHOD -> ENSURE SOMETHING ADDED TO EVERY NETWORK, and EVERY SEED OPTION
    def essential(self):
        #TURN OFF AND ON SPLIT HALF ANALYSIS
        #TURN OFF AND ON GENETIC ALGORITHM (SPECIFY PARAMETERS)
        #SPECIFY NUMBER OF PIPELINES
        #DATA INPUT
        #CACHE INTERMEDIATES
        #SETUP SEEDS
        #LINEAR, MULTI, or CLUSTER RUN
        #RUN ANALYSIS OR SAVE CONFIG FILE
        #HASHING TYPE #BLOCK GENERATE DICT IF VALUES NOT INPUTTED
        self.aesthetic_frame_main = ttk.Frame(self.tab_main)
        self.aesthetic_frame_main.grid(row=0, column=0, sticky="nsew")
        self.aesthetic_frame_main.columnconfigure(0, weight=1)
        ttk.Label(self.aesthetic_frame_main, text='Essential Configurations', font='Helvetica 12 bold').grid(row=0)
        self.main_config = ttk.Frame(self.aesthetic_frame_main, relief='groove', padding='0.5i')
        self.main_config.grid(row=1)
        
        #Number of networks to target
        #Types of seeding to try
        ttk.Label(self.main_config, text='Seed types:', font='Helvetica 12').grid(row=0)
        self.atlas = tk.IntVar()
        self.atlas_check = tk.Checkbutton(self.main_config, text='Harvard-Oxford Atlas/Saved Mask', variable=self.atlas)
        self.atlas_check.grid(row=2, column=0)
        self.data = tk.IntVar()
        self.data_check = tk.Checkbutton(self.main_config, text='ReHo Defined Seed', variable=self.data)
        self.data_check.grid(row=2, column=1)
        self.coords = tk.IntVar()
        self.coords_check = tk.Checkbutton(self.main_config, text='User Defined ROI', variable=self.coords)
        self.coords_check.grid(row=2, column=2)
        
        ttk.Label(self.main_config, text='Number of networks:', font='Helvetica 12').grid(row=3)
        self.networks = tk.Entry(self.main_config)
        self.networks.insert(4, "1")
        self.networks.grid(row=3, column=1)
        
        self.seedbutton = tk.Button(self.main_config, text='Define Seeds', command=self.define_seeds)
        self.seedbutton.grid(row=4, column=0)
        
    
    def define_seeds(self):#for adding a new one: type in box, hit enter, add to dictionary, clear, for new file, just select multiple files (post processing is copying them to local location)
        self.define_frame = ttk.Frame(self.tab_main)
        self.define_frame.grid(row=0, column=0, sticky="nsew")
        self.define_frame.columnconfigure(0, weight=1)
        ttk.Label(self.define_frame, text='Seed Configuration', font='Helvetica 12 bold').grid(row=0)
        
        try:
            networks = int(self.networks.get())
        except:
            raise ValueError('Invalid number of networks. Please enter an integer - defaulting to 1 network')
            networks = 1
            
        #labels = ['Brain Region (can use abbreviation but must']
        modes = [self.atlas.get(), self.data.get(), self.coords.get()]
        
        if not sum(modes):
            raise ValueError('At least one seed definition type must be selected')
        
        counter = 1
        self.out_dic['level1']['~construct~Finfo_rest_type']['gene'] = []
        for val in range(networks):
            ttk.Label(self.define_frame, text='Network {val}:'.format(val=val), font='Helvetica 12 bold').grid(row=counter+val)
            vars(self)['define_network'+str(val)] = ttk.Frame(self.define_frame, relief='groove', padding='0.5i')
            vars(self)['define_network'+str(val)].grid(row=counter+val+1)
            for i, mode in enumerate(modes):
                if mode:
                    if i == 0:
                        vars(self)['atlas_dropdown'+str(val)] = []
                        fsl = os.getenv('FSLDIR')
                        for atlas in glob.glob(fsl + '/data/atlases/HarvardOxford*.xml'):
                            tree = ET.parse(atlas)
                            root = tree.getroot()
                        
                            for label in root.iter('label'):
                                vars(self)['atlas_dropdown'+str(val)].append(label.text)
                        
                        self.out_dic['level1']['~construct~Finfo_rest_type']['gene'].append(i)
                        
                        ttk.Label(vars(self)['define_network'+str(val)], text='Harvard-Oxford Atlas Brain Region:').grid(row=counter+val)
                        vars(self)['drop_select'+str(val)] = tk.StringVar()
                        vars(self)['drop_select'+str(val)].set("")
                        vars(self)['option_menu'+str(val)] = tk.OptionMenu(vars(self)['define_network'+str(val)], vars(self)['drop_select'+str(val)], *vars(self)['atlas_dropdown'+str(val)])
                        vars(self)['option_menu'+str(val)].grid(row=counter+val, column=1)
                        counter += 1
                        
                        ttk.Label(vars(self)['define_network'+str(val)], text='Minimum threshold:').grid(row=counter+val)
                        vars(self)['atlas_min'+str(val)] = tk.Entry(vars(self)['define_network'+str(val)])
                        vars(self)['atlas_min'+str(val)].insert(4, '0')
                        vars(self)['atlas_min'+str(val)].grid(row=counter+val, column=1)
                        counter += 1
                        
                        ttk.Label(vars(self)['define_network'+str(val)], text='Maximum threshold').grid(row=counter+val)
                        vars(self)['atlas_max'+str(val)] = tk.Entry(vars(self)['define_network'+str(val)])
                        vars(self)['atlas_max'+str(val)].insert(4, '95')
                        vars(self)['atlas_max'+str(val)].grid(row=counter+val, column=1)
                        counter += 1
                        
                        vars(self)['called_'+str(i)+'_'+str(val)] = 0
                        ttk.Label(vars(self)['define_network'+str(val)], text='Multiverse Atlas ROIs:').grid(row=counter+val)
                        vars(self)['called_'+str(i)+'_'+str(val)+'label'] = ttk.Label(vars(self)['define_network'+str(val)], text=str(vars(self)['called_'+str(i)+'_'+str(val)]))
                        vars(self)['called_'+str(i)+'_'+str(val)+'label'].grid(row=counter+val, column=1)
                        counter += 1
                        
                        vars(self)['atlas_add'+str(val)] = tk.Button(vars(self)['define_network'+str(val)], text='Add Brain Region', command=partial(self.add_option, i, val, 'called_'+str(i)+'_'+str(val), True, networks))
                        vars(self)['atlas_add'+str(val)].grid(row=counter+val, column=1)
                        #counter += 1
                        
                        vars(self)['atlas_add'+str(val)] = tk.Button(vars(self)['define_network'+str(val)], text='Add Saved Mask', command=partial(self.add_option, i, val, 'called_'+str(i)+'_'+str(val), False, networks))
                        vars(self)['atlas_add'+str(val)].grid(row=counter+val, column=0)
                        counter += 1
        
        self.define_frame.tkraise()
        
    def add_option(self, i, val, called_name, atlas, networks):
        if i == 0:
            options = len(self.out_dic['level1']['~construct~Finfo_rest_seedinfo']) - 1
            if atlas:
                mini = vars(self)['atlas_min'+str(val)].get()
                if re.search('[^a-zA-Z0-9\.]+', mini) or not mini:
                    raise ValueError('Expected an integer, but {mini} was provided instead. Using 0'.format(mini=mini))
                    mini = 0
                else:
                    mini = self.to_num(mini)
                    if mini > 95:
                        raise ValueError('Value >95 provided, which will result in zero map. Using 95')
                        mini = 95
                    elif mini < 0:
                        raise ValueError('Value <0 provided, using 0')
                        mini = 0
                        
                        
                maxi = vars(self)['atlas_max'+str(val)].get()
                if re.search('[^a-zA-Z0-9\.]+', maxi) or not maxi:
                    raise ValueError('Expected an integer, but {maxi} was provided instead. Using 95'.format(maxi=maxi))
                    maxi = 95
                else:
                    maxi = self.to_num(maxi)
                    if maxi > 100:
                        raise ValueError('Value >100 provided, which will result in zero map. Using 99')
                        maxi = 95
                    elif maxi < 0:
                        raise ValueError('Value <0 provided, using 0')
                        maxi = 0
                        
                maxi = self.to_num(vars(self)['atlas_max'+str(val)].get())
                if maxi < mini:
                    seedinfo = (vars(self)['drop_select'+str(val)].get(), mini)
                else:
                    seedinfo = (vars(self)['drop_select'+str(val)].get(), mini, maxi)
            else:
                seedinfo = [tk.filedialog.askopenfilename()] #if filenames list() instead of []
                dir = os.path.dirname(__file__)
                for l, seed in enumerate(seedinfo):
                    file_name = re.search('.*/(.*)', seed).group(1)
                    destination = os.path.join(dir, '..', 'seed_masks')
                    if not os.path.isdir(destination):
                        os.mkdir(destination)
                    
                    copy2(seed, os.path.join(destination, file_name))
                    seedinfo[l] = os.path.join('..', 'seed_masks', file_name)
                    
                seedinfo = seedinfo[0]
                
                
            if vars(self)[called_name] < options:
                for j in range(vars(self)[called_name], options):#TEST THIS-> WAS OPTIONS BUT WAS BROKEN
                    if len(self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][j]) <= val:
                        self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][j].append(seedinfo)
                    else:
                        self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][val] = seedinfo
            else:
                if options > 0:
                    self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options] = self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options-1].copy()
                    self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options][val] = seedinfo
                else:
                    self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options] = [seedinfo] * networks
        
            self.out_dic['level1']['~construct~Finfo_rest_seedinfo']['genes'] = {'low': -0.5, 'high': options-0.5}
            vars(self)['drop_select'+str(val)].set(" ")
            vars(self)['atlas_max'+str(val)].delete(0, 'end')
            vars(self)['atlas_max'+str(val)].insert(4, '95')
            vars(self)['atlas_min'+str(val)].delete(0, 'end')
            vars(self)['atlas_min'+str(val)].insert(4, '0')
            vars(self)[called_name] += 1
            vars(self)[called_name+'label']['text'] = vars(self)[called_name]
    
    def get_format(self):
        with open('/Users/grahamseasons/fMRI/analysis_multiverse/default.json') as f:
            data = json.load(f)
        
        return data
    
    def strings(self, stage, alias, default, value_map, param_name_full, value, butt_text, counter):
        vars(self)['aesthetic_frame_'+param_name_full] = ttk.Frame(vars(self)['tab_window_'+stage], relief="groove", padding="0.5i")
        vars(self)['aesthetic_frame_'+param_name_full].grid(row=counter, column=0, sticky="nsew")
        vars(self)['aesthetic_frame_'+param_name_full].columnconfigure(0, weight=1)
        
        ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text=alias).grid(row=1)
        #CONSTRUCT NODE PARAMETER NAME AND ADD TO DICTIONARY
        ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text='Value Pool:').grid(row=2)
        vars(self)['valuepool_'+param_name_full] = ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text=str(value_map))
        vars(self)['valuepool_'+param_name_full].grid(row=2, column=1)

        vars(self)[param_name_full+'configbutton'] = tk.Button(vars(self)['aesthetic_frame_'+param_name_full], text='Configure', 
                                                                command=partial(self.mapped_config_frame, stage, param_name_full, alias, value, value_map))
        vars(self)[param_name_full+'configbutton'].grid(row=3, column=1)
        
        vars(self)[param_name_full+'onbutton'] = tk.Button(vars(self)['aesthetic_frame_'+param_name_full], text=butt_text, command=partial(self.change, param_name_full+'onbutton', 'valuepool_'+param_name_full, default, value, value_map, stage, param_name_full, param_name_full+'configbutton'))
        vars(self)[param_name_full+'onbutton'].grid(row=1, column=1)
        
        
    
    def mapped_config_frame(self, stage, param_name_full, alias, value, value_map):
        vars(self)[param_name_full+'frame'] = tk.Frame(vars(self)['tab_'+stage])
        vars(self)[param_name_full+'frame'].grid(row=0, column=0, sticky="nsew")
        vars(self)[param_name_full+'frame'].columnconfigure(0, weight=1)
        ttk.Label(vars(self)[param_name_full+'frame'], text=alias + ' Configuration', font='Helvetica 12 bold').grid(row=0)
        
        for val in value:
            vars(self)[param_name_full+str(val)+'_var'] = tk.IntVar()
            vars(self)[param_name_full+str(val)] = tk.Checkbutton(vars(self)[param_name_full+'frame'], text=str(value_map[val]), variable=vars(self)[param_name_full+str(val)+'_var'])
            vars(self)[param_name_full+str(val)].grid(row=val+1)
            
        vars(self)[param_name_full+'frame_button'] = tk.Button(vars(self)[param_name_full+'frame'], text='Exit', 
                                                                   command=partial(self.ret, stage))
        vars(self)[param_name_full+'frame_button'].grid(row=val+2, column=0)
        
        vars(self)[param_name_full+'frame_button'] = tk.Button(vars(self)[param_name_full+'frame'], text='Set', 
                                                               command=partial(self.selection, value, value_map, param_name_full, stage))
        vars(self)[param_name_full+'frame_button'].grid(row=val+2, column=1)
        
        vars(self)[param_name_full+'frame'].tkraise()
        
            
    def selection(self, value, value_map, param_name_full, stage):
        lst = []
        updated = []
        for val in value:
            temp = vars(self)[param_name_full+str(val)+'_var'].get()
            if temp:
                lst.append(val)
                updated.append(value_map[val])
        if not lst:
            raise ValueError('No values selected. Must select 1 or more')
            
        self.out_dic[stage][param_name_full]['gene'] = lst
        vars(self)['valuepool_'+param_name_full]['text'] = str(updated)
        
        vars(self)['tab_window_'+stage].tkraise()
        
    
    def numerical(self, stage, alias, default, value_map, param_name_full, value, butt_text, counter):
        vars(self)['aesthetic_frame_'+param_name_full] = ttk.Frame(vars(self)['tab_window_'+stage], relief="groove", padding="0.5i")
        vars(self)['aesthetic_frame_'+param_name_full].grid(row=counter, column=0, sticky="nsew")
        vars(self)['aesthetic_frame_'+param_name_full].columnconfigure(0, weight=1)
        
        ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text=alias).grid(row=1)
        #CONSTRUCT NODE PARAMETER NAME AND ADD TO DICTIONARY
        ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text='Value Pool:').grid(row=2)
        vars(self)['valuepool_'+param_name_full] = ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text=str(value))
        vars(self)['valuepool_'+param_name_full].grid(row=2, column=1)

        ttk.Label(vars(self)['aesthetic_frame_'+param_name_full], text='Input Type:').grid(row=3)
        vars(self)[param_name_full+'rangebutton'] = tk.Button(vars(self)['aesthetic_frame_'+param_name_full], text='Range', command=partial(self.change, param_name_full+'rangebutton', 'valuepool_'+param_name_full, default, value, value_map, stage, param_name_full, ''))
        vars(self)[param_name_full+'rangebutton'].grid(row=3, column=1)
        
        vars(self)[param_name_full+'configbutton'] = tk.Button(vars(self)['aesthetic_frame_'+param_name_full], text='Configure', 
                                                                command=partial(self.config_frame, param_name_full+'rangebutton', stage, param_name_full, alias))
        vars(self)[param_name_full+'configbutton'].grid(row=4, column=1)
        
        vars(self)[param_name_full+'onbutton'] = tk.Button(vars(self)['aesthetic_frame_'+param_name_full], text=butt_text, command=partial(self.change, param_name_full+'onbutton', 'valuepool_'+param_name_full, default, value, value_map, stage, param_name_full, param_name_full+'configbutton'))
        vars(self)[param_name_full+'onbutton'].grid(row=1, column=1)

        
    #OFF ALWAYS GOES TO MY DEFAULT -> ON GO TO USER SET VALUE
    #ON/OFF HERE NEEDS TO REFER TO THE DICTIONARY, NOT SOLELY JSON FILE
    def change(self, name, def_name, default, value, value_map, stage, param_name_full, configure):
        if vars(self)[name]['text'] == 'On':
            vars(self)[name].config(text='Off')
            vars(self)[configure]['state'] = 'disabled'
            if value_map:
                vars(self)[def_name].config(text=str(value_map[default]))
                self.out_dic[stage][param_name_full]['gene'] = default
            else:
                vars(self)[def_name].config(text=str(value[default]))
                self.out_dic[stage][param_name_full]['gene'] = value[default]
        elif vars(self)[name]['text'] == 'Off':
            vars(self)[name].config(text='On')
            vars(self)[configure]['state'] = 'normal'
            if value_map:
                vars(self)[def_name].config(text=str(value_map))
                self.out_dic[stage][param_name_full]['gene'] = value #THIS NEEDS TO BE IN THE FORM OF NUMBERS, NOT STRINGS
            else:
                vars(self)[def_name].config(text=str(value))
                self.out_dic[stage][param_name_full]['gene'] = value
        elif vars(self)[name]['text'] == 'Range':
            vars(self)[name].config(text='Manual')
        elif vars(self)[name]['text'] == 'Manual':
            vars(self)[name].config(text='Range')
                
    def config_frame(self, range_state, stage, param_name_full, alias):
        vars(self)[param_name_full+'frame'] = tk.Frame(vars(self)['tab_'+stage])
        vars(self)[param_name_full+'frame'].grid(row=0, column=0, sticky="nsew")
        vars(self)[param_name_full+'frame'].columnconfigure(0, weight=1)
        ttk.Label(vars(self)[param_name_full+'frame'], text=alias + ' Configuration', font='Helvetica 12 bold').grid(row=0)
        if vars(self)[range_state]['text'] == 'Range':
            ttk.Label(vars(self)[param_name_full+'frame'], text='Low: ').grid(row=1)
            vars(self)[param_name_full+'low'] = tk.Entry(vars(self)[param_name_full+'frame'])
            vars(self)[param_name_full+'low'].grid(row=1, column=1)
            ttk.Label(vars(self)[param_name_full+'frame'], text='High: ').grid(row=2)
            vars(self)[param_name_full+'high'] = tk.Entry(vars(self)[param_name_full+'frame'])
            vars(self)[param_name_full+'high'].grid(row=2, column=1)
            ttk.Label(vars(self)[param_name_full+'frame'], text='Step: ').grid(row=3)
            vars(self)[param_name_full+'step'] = tk.Entry(vars(self)[param_name_full+'frame'])
            vars(self)[param_name_full+'step'].grid(row=3, column=1)
            
            vars(self)[param_name_full+'frame_button'] = tk.Button(vars(self)[param_name_full+'frame'], text='Exit', 
                                                                   command=partial(self.ret, stage))
            vars(self)[param_name_full+'frame_button'].grid(row=4, column=0)
            
            vars(self)[param_name_full+'frame_button'] = tk.Button(vars(self)[param_name_full+'frame'], text='Set', 
                                                                   command=partial(self.set_values, range_state, param_name_full, stage))
            vars(self)[param_name_full+'frame_button'].grid(row=4, column=1)
        elif vars(self)[range_state]['text'] == 'Manual':
            ttk.Label(vars(self)[param_name_full+'frame'], text='Values (comma separated list): ').grid(row=1)
            vars(self)[param_name_full+'manual'] = tk.Entry(vars(self)[param_name_full+'frame'])
            vars(self)[param_name_full+'manual'].grid(row=2)
            
            vars(self)[param_name_full+'frame_button'] = tk.Button(vars(self)[param_name_full+'frame'], text='Exit', 
                                                                   command=partial(self.ret, stage))
            vars(self)[param_name_full+'frame_button'].grid(row=3, column=0)
            
            vars(self)[param_name_full+'frame_button'] = tk.Button(vars(self)[param_name_full+'frame'], text='Set', 
                                                                   command=partial(self.set_values, range_state, param_name_full, stage))
            vars(self)[param_name_full+'frame_button'].grid(row=3, column=1)
            
        vars(self)[param_name_full+'frame'].tkraise()
            
    def set_values(self, range_state, param_name_full, stage):
        if vars(self)[range_state]['text'] == 'Range':
            low = self.to_num(vars(self)[param_name_full+'low'].get())
            high = self.to_num(vars(self)[param_name_full+'high'].get())
            step = self.to_num(vars(self)[param_name_full+'step'].get())
            if step:
                self.out_dic[stage][param_name_full]['gene'] = {'low': low, 'high': high, 'step': step}
                updated = {'low': low, 'high': high, 'step': step}
            else:
                self.out_dic[stage][param_name_full]['gene'] = {'low': low, 'high': high}
                updated = {'low': low, 'high': high}
            
            vars(self)['valuepool_'+param_name_full]['text'] = str(step)
        elif vars(self)[range_state]['text'] == 'Manual':
            updated = vars(self)[param_name_full+'manual'].get()
            updated = updated.replace(' ', '').split(',')
            for i, val in enumerate(updated):
                updated[i] = self.to_num(val)
            
            self.out_dic[stage][param_name_full]['gene'] = updated
            
        vars(self)['valuepool_'+param_name_full]['text'] = str(updated)
        vars(self)['tab_window_'+stage].tkraise()
    
    def ret(self, stage):
        vars(self)['tab_window_'+stage].tkraise()
            
    def to_num(self, str_):
        try:
            return int(str_)
        except:
            try:
                return float(str_)
            except:
                raise ValueError('Invalid input type, must be int or float')
                
    def get_full_mult(self, stage, node, key, toggle, orig_entry):
        A = 3
        


MultiverseConfig()


start = argparse.ArgumentParser()
start.add_argument('-c', '--config', action='store_true', help='configure multiverse parameter file')
start.add_argument('-g', '--gui', action='store_true', help='run graphical user interface')
start.add_argument('-r', '--run', action='store_true', help='run multiverse analysis')
start.add_argument('-d', '--data', type=str, metavar="DATA_DIR", action='store', help='path to BIDS formatted data directory')
start.add_argument('-o', '--out', type=str, metavar="OUT_DIR", action='store', help='path to store outputs')
#SLURM, MULTIPROC
args = start.parse_args()

if args.data == None and not args.gui:
    print('Either -d must be specified, or -g must be used')
    sys.exit()
elif not os.path.isdir(args.data):
    print('The specified data path does not exist')
    sys.exit()
    
if args.out == None and not args.gui:
    print('Either -o must be specified, or -g must be used')
    sys.exit()
elif not os.path.isdir(args.out):
    print('The specified data path does not exist')
    print('Creating directory at specified path')
    os.makedirs(args.out)
    
if args.out == None and not args.gui:
    print('Either -o must be specified, or -g must be used')
    sys.exit()
elif not os.path.isdir(args.out):
    print('The specified data path does not exist')
    print('Creating directory at specified path')
    os.makedirs(args.out)

if not args.config:
    USEDEFAULT=True #LOOK FOR FILE IN DEFAULT LOCATION, OR GENERATE
    #STILL NEED TO ASK FOR COORDINATES AND MASKS
elif args.gui:
    USEGUI=True
else:
    USECLI=True
    

#HAS TO BE GUI FOR SETUP, OR MANUAL
# =============================================================================
# index = 0
# while True:
#     gene = genes[index]
#     name = list(gene.keys())[0]
#     value = list(gene.values())
#     print("Parameter: {name}\n".format(name=name))
#     print("Default range: {value}\n".format(value=value))
#     while True:
#         switch = input("Include in multiverse analysis [y/n]? ")
#         if switch == 'y':
#             switch = True
#             newprompt = "Input new range (keep same format as default): "
#             break
#         elif switch == 'n':
#             newprompt = "Input constant value: "
#             switch = False
#             break
#         elif switch == 'exit':
#             sys.exit()
#         else:
#             switch = input("Invalid input {switch}! Please try again, specifying either 'y' or 'n' ")
#     
#     if switch:
#         
# =============================================================================

#Probably need to include links, but have flag to not display in gui
create_docker_command = '-v {code}:/root/multiverse/code -v {data}:/data -v {work_dir}:/scratch -w scratch'
create_docker_command = create_docker_command.format(code=code, data=data_dir, work_dir=exp_dir)
#THIS SHOULD SET UP DOCKER, AND ANALYSIS CODE, DOCKER SHOULD THEN AUTOMATICALLY RUN ANALYSIS CODE
#SAVE CONFIG FILE WITH DEFAULTS AND ALTERED VALUES TO LOCATION
subprocess.run([])

# =============================================================================
# neurodocker generate docker --base debian:stretch --pkg-manager apt --fsl version=6.0.0 --ants version=2.3.1 --afni version=latest --user multiverse --miniconda miniconda_version="4.3.31" create_env=multiverse conda_install="python=3.8 pytest traits pandas matplotlib scikit-learn scikit-image seaborn numpy" pip_install="nipype nilearn nibabel niworkflows pygad kivy" --user root --install nano --volume /data /root/multiverse/code /scratch -w /scratch > multiverse.Dockerfile
# neurodocker generate docker --base debian:stretch --pkg-manager apt --fsl version=5.0.11 --ants version=2.3.1 --afni version=latest --miniconda miniconda_version="4.3.31" create_env=multiverse conda_install="python=3.8 pytest traits pandas matplotlib scikit-learn scikit-image seaborn numpy" pip_install="nipype nilearn nibabel niworkflows pygad kivy" --user root --install nano --user multiverse --volume /data /root/multiverse/code /scratch -w /scratch > /Volumes/NewVolume/docker/multiverse.Dockerfile
# docker build --tag multiverse --file multiverse.Dockerfile .
# =============================================================================
#docker run -v /Users/grahamseasons/fMRI/analysis_multiverse/updated:/root/multiverse/code -v /Volumes/NewVolume/super_agers:/data -v /Volumes/NewVolume/docker_test:/scratch -w /scratch multiverse