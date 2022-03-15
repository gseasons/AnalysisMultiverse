#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 14:48:48 2021

@author: grahamseasons
"""
import os, re, glob
from shutil import copy2
import xml.etree.ElementTree as ET
import json, pickle
import tkinter as tk
from tkinter import ttk
from functools import partial
import numpy as np

class AutoScrollbar(tk.Scrollbar):
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        tk.Scrollbar.set(self, lo, hi)
    def pack(self, **kw):
        raise tk.TclError("cannot use pack with this widget")
    def place(self, **kw):
        raise tk.TclError("cannot use place with this widget")

class MultiverseConfig():
    def __init__(self, rerun, data_dir="N/A", output="N/A"):
        self.data_dir = data_dir
        self.output = output
        self.master = tk.Tk()
        self.master.geometry("900x550")
        
        self.main_set_butt = {'state': 'deactivate'}
        
        self.vscroll = AutoScrollbar(self.master)
        self.vscroll.grid(row=0, column=1, sticky='ns')
        self.canvas = tk.Canvas(self.master, yscrollcommand=self.vscroll.set)
        self.canvas.grid(row=0, column=0, sticky='nsew')
        self.canvas.rowconfigure(0, weight=1)
        self.vscroll.config(command=self.canvas.yview)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)
        
        self.out_dic = {}
        self.configure = {'split_half': False, 'debug': False, 'rerun': rerun}
        
        self.tabcontrol = ttk.Notebook(self.canvas)
        self.tab_main = ttk.Frame(self.tabcontrol)
        self.tab_main.columnconfigure(0, weight=1)
        self.tabcontrol.add(self.tab_main, text='Essential')
        
        label_dic = {'preprocess': 'Preprocessing', 'level1': 'Level 1 Analysis', 'level2': 'Level 2 Analysis', 'level3': 'Level 3 Analysis', 'correction': 'Multiple Comparisons Correction'}
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
                if node['alias'] != 'link':
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
                        
                    if on_off == "on":
                        butt_text = 'On'
                        self.out_dic[stage][param_name_full] = {'gene': value}
                    else:
                        butt_text = 'Off'
                        self.out_dic[stage][param_name_full] = {'gene': [value[default]]}
                        
                    if 'show' in param:
                        if value_map:
                            for i in value:
                                self.out_dic[stage][param_name_full][i] = value_map[i]
                        continue
                        
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
            
        self.essential()
        self.tabcontrol.pack(expand=1, fill="both")
        self.canvas.iid = self.canvas.create_window(0, 0, anchor='nw', window=self.tabcontrol)
        self.canvas.bind("<Configure>", self.canvas_configure)
        self.canvas.bind_all('<MouseWheel>', lambda e: self.scroll(e))
        
        self.tabcontrol.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox('all'))
        self.master.mainloop()
        
    def canvas_configure(self, event):
        canvas = event.widget
        canvas.itemconfigure(canvas.iid, width=canvas.winfo_width())
        
    def scroll(self, event):
        self.canvas.yview_scroll(-1*event.delta, tk.UNITS)
    
    def essential(self):
        self.aesthetic_frame_main = ttk.Frame(self.tab_main)
        self.aesthetic_frame_main.grid(row=0, column=0, sticky="nsew")
        self.aesthetic_frame_main.columnconfigure(0, weight=1)
        ttk.Label(self.aesthetic_frame_main, text='Essential Configurations', font='Helvetica 12 bold').grid(row=0)
        self.main_config = ttk.Frame(self.aesthetic_frame_main, relief='groove', padding='0.5i')
        self.main_config.grid(row=1)

        ttk.Label(self.main_config, text='Seed types:', font='Helvetica 12').grid(row=0)
        self.atlas = tk.IntVar()
        self.atlas_check = tk.Checkbutton(self.main_config, text='Harvard-Oxford Atlas/Saved Mask', variable=self.atlas, command=self.allow)
        self.atlas_check.grid(row=2, column=0)
        self.data = tk.IntVar()
        self.data_check = tk.Checkbutton(self.main_config, text='ReHo Defined Seed', variable=self.data, command=self.allow)
        self.data_check.grid(row=2, column=1)
        self.coords = tk.IntVar()
        self.coords_check = tk.Checkbutton(self.main_config, text='User Defined ROI', variable=self.coords, command=self.allow)
        self.coords_check.grid(row=2, column=2)
        
        ttk.Label(self.main_config, text='Number of networks:', font='Helvetica 12').grid(row=3)
        self.networks = tk.Entry(self.main_config)
        self.networks.insert(4, "1")
        self.networks.grid(row=3, column=1)
        self.networks.bind("<Key>", self.valid_networks)
        
        self.seedbutton = tk.Button(self.main_config, text='Define Seeds', command=self.define_seeds)
        self.seedbutton['state'] = 'disabled'
        self.seedbutton.grid(row=4, column=0)
        
        ttk.Label(self.main_config, text='Number of pipelines:', font='Helvetica 12').grid(row=5)
        self.pipelines = tk.Entry(self.main_config)
        self.pipelines.insert(4, "50")
        self.pipelines.grid(row=5, column=1)
        self.pipelines.bind("<Key>", self.valid_pipelines)
        
        ttk.Label(self.main_config, text='Genetic Algorithm:', font='Helvetica 12').grid(row=6)
        self.ga = tk.Button(self.main_config, text='Off', command=self.genetic_switch)
        self.ga.grid(row=6, column=1)
        
        self.ga_config = tk.Button(self.main_config, text='Genetic Algorithm Config', command=self.genetic_config)
        self.ga_config['state'] = 'disabled'
        self.ga_config.grid(row=7, column=0)
        
        ttk.Label(self.main_config, text='Split-Half Analysis:').grid(row=8)
        self.split = ttk.Label(self.main_config, text='Off')
        self.split.grid(row=8, column=1)
        
        ttk.Label(self.main_config, text='Debug mode (cached intermediates):').grid(row=9)
        self.debug_butt = ttk.Button(self.main_config, text='Off', command=self.debug)
        self.debug_butt.grid(row=9, column=1)
        
        ttk.Label(self.main_config, text='Processing mode (SLURM for compute cluster):').grid(row=10)
        self.process_mode_var = tk.StringVar()
        self.process_mode_var.set('MultiProc')
        self.parent_selection = tk.OptionMenu(self.main_config, self.process_mode_var, 'MultiProc', 'Linear', 'SLURM', command=self.slurm_block)
        self.parent_selection.grid(row=10, column=1)
        
        ttk.Label(self.main_config, text='Maximum Output Storage (GB, ignored for debug mode):', font='Helvetica 12').grid(row=11)
        self.storage = tk.Entry(self.main_config)
        self.storage.insert(4, "20000")
        self.storage.grid(row=11, column=1)
        
        ttk.Label(self.main_config, text='Data directory:').grid(row=12)
        ttk.Label(self.main_config, text=self.data_dir).grid(row=12, column=1)
        ttk.Label(self.main_config, text='Output directory:').grid(row=13)
        ttk.Label(self.main_config, text=self.output).grid(row=13, column=1)
        
        self.return_buttons = ttk.Frame(self.aesthetic_frame_main, relief='groove', padding='0.5i')
        self.return_buttons.grid(row=2)
        
        self.configure_now = tk.Button(self.return_buttons, text='Save Configuration', command=partial(self.done, False))
        self.configure_now['state'] = 'disabled'
        self.configure_now.grid(row=0, column=0)
        
        self.run_now = tk.Button(self.return_buttons, text='Run', command=partial(self.done, True))
        self.run_now['state'] = 'disabled'
        self.run_now.grid(row=0, column=1)
        
    def valid_networks(self, change):
        try:
            if change.char == '\x7f':
                new_pipe = int(self.networks.get()[:-1])
            else:
                new_pipe = int(self.networks.get() + change.char)
                
            if new_pipe < 1:
                raise ValueError('Number of networks must be a positive integer')
                
            self.seedbutton['state'] = 'normal'
        except ValueError:
            self.seedbutton['state'] = 'disabled'
            
    def valid_pipelines(self, change):
        try:
            if change.char == '\x7f':
                new_pipe = int(self.pipelines.get()[:-1])
            else:
                new_pipe = int(self.pipelines.get() + change.char)
                
            if new_pipe < 1:
                raise ValueError('Number of pipelines must be a positive integer')
                
            self.ga['state'] = 'normal'
            if self.ga['text'] == 'On':
                self.ga_config['state'] = 'normal'
            
            self.validate()
        except ValueError:
            self.configure_now['state'] = 'disabled'
            self.run_now['state'] = 'disabled'
            self.ga_config['state'] = 'disabled'
            self.ga['state'] = 'disabled'
        
    def validate(self):
        if 'networks' in self.configure:
            if 'num_generations' in self.configure or self.ga['text'] == 'Off':
                self.configure_now['state'] = 'normal'
                if self.data_dir != 'N/A' and self.output != 'N/A':
                    self.run_now['state'] = 'normal'
        
    def slurm_block(self, change):
        if self.process_mode_var.get() == 'SLURM':
            self.run_now['state'] = 'disabled'
            self.slurm_config()
        elif self.main_set_butt['state'] == 'normal':
            if 'num_generations' in self.configure or self.ga['text'] == 'Off':
                if self.data_dir != 'N/A' and self.output != 'N/A':
                    self.run_now['state'] = 'normal'
            else:
                self.run_now['state'] = 'disabled'
                
    def slurm_config(self):
        self.slurm_frame_main = ttk.Frame(self.tab_main)
        self.slurm_frame_main.grid(row=0, column=0, sticky="nsew")
        self.slurm_frame_main.columnconfigure(0, weight=1)
        ttk.Label(self.slurm_frame_main, text='Batch Submission Settings', font='Helvetica 12 bold').grid(row=0)
        
        self.slurm_frame = ttk.Frame(self.slurm_frame_main, relief='groove', padding='0.5i')
        self.slurm_frame.grid(row=1, column=0, sticky="nsew")
        self.slurm_frame.columnconfigure(0, weight=1)
        
        ttk.Label(self.slurm_frame, text='SLURM account:').grid(row=1)
        self.account = ttk.Entry(self.slurm_frame)
        self.account.insert(4, 'def-')
        self.account.grid(row=1, column=1)
        
        ttk.Label(self.slurm_frame, text='Number of CPUs:').grid(row=2)
        self.cpus = ttk.Entry(self.slurm_frame)
        self.cpus.insert(4, self.pipelines.get())
        self.cpus.grid(row=2, column=1)
        
        ttk.Label(self.slurm_frame, text='CPUs per node (on compute cluster):').grid(row=3)
        self.nodes = ttk.Entry(self.slurm_frame)
        self.nodes.insert(4, '32')
        self.nodes.grid(row=3, column=1)
        
        ttk.Label(self.slurm_frame, text='Pipelines per batch:').grid(row=4)
        self.batches = ttk.Entry(self.slurm_frame)
        self.batches.insert(4, str(np.ceil(int(self.pipelines.get())/4).astype(int)))
        self.batches.grid(row=4, column=1)
        
        ttk.Label(self.slurm_frame, text='Memory required per CPU (MB):').grid(row=5)
        self.mem = ttk.Entry(self.slurm_frame)
        self.mem.insert(4, '6000')
        self.mem.grid(row=5, column=1)
        #NEED BETTER BENCHMARK
        ttk.Label(self.slurm_frame, text='Time ~2H * (subject,scan,pipeline) / CPUs').grid(row=6)
        
        ttk.Label(self.slurm_frame, text='Job time (Days-Hours:Mins):').grid(row=7)
        self.time = ttk.Entry(self.slurm_frame)
        self.time.insert(4, '0-00:00')
        self.time.grid(row=7, column=1)
        
        self.set_genetic = tk.Button(self.slurm_frame, text='Set', command=partial(self.set_slurm, False))
        self.set_genetic.grid(row=8)
        
        self.slurm_frame_main.tkraise()
        
    def set_slurm(self, auto):
        node_request = int(self.cpus.get()) / int(self.nodes.get())
        self.configure['nodes'] = str(np.ceil(node_request).astype(int))
        self.configure['ntasks'] = str(np.ceil(node_request).astype(int) * int(self.nodes.get()))
        self.configure['cpu_node'] = self.nodes.get()
        self.configure['batches'] = int(self.batches.get())
        self.configure['account'] = self.account.get()
        self.configure['mem'] = self.mem.get()
        self.configure['time'] = self.time.get()
        
        self.aesthetic_frame_main.tkraise()
        
    def done(self, run):
        self.configure['processing'] = self.process_mode_var.get()
        self.configure['storage'] = float(self.storage.get())
        if 'pipelines' not in self.configure:
            self.configure['pipelines'] = int(self.pipelines.get())
        
        dir = os.path.dirname(os.path.abspath(__file__))
        settings = 'multiverse_configuration.pkl'
        general = 'general_configuration.pkl'
        
        destination = os.path.join(dir, '../configuration')
        
        if not os.path.isdir(destination):
            os.mkdir(destination)
            
        settings = os.path.join(destination, settings)
        general = os.path.join(destination, general)
        reformatted = self.format_out(self.out_dic)
        self.save(settings, reformatted)
        self.save(general, self.configure)
        #RESUME MODE
        self.run_now = run
        
        self.master.destroy()
        self.master.quit()
        
    def format_out(self, out):
        new = []
        prior = ''
        for level in out:
            if prior:
                new.append({'end_'+prior: level})
            for key in out[level]:
                values = out[level][key]
                if values['gene']:
                    if len(values) > 1:
                        temp = {key: values['gene']}
                        temp.update(values)
                        temp.pop('gene')
                        new.append(temp)
                    else:
                        new.append({key: values['gene']})
            
            prior = level
            
        new.append({'end_'+prior: 'end'})
        
        return new
            
    def save(self, file_name, dictionary):
        file = open(file_name, 'wb')
        pickle.dump(dictionary, file)
        file.close()
        
    def debug(self):
        if self.debug_butt['text'] == 'On':
            self.debug_butt['text'] = 'Off'
            self.configure['debug'] = False
        else:
            self.debug_butt['text'] = 'On'
            self.configure['debug'] = True
            
    def allow(self):
        if int(self.networks.get()) > 0 and (self.data.get() or self.coords.get() or self.atlas.get()):
            self.seedbutton['state'] = 'normal'
        else:
            self.seedbutton['state'] = 'disabled'
    
    def genetic_switch(self):
        if self.ga['text'] == 'Off':
            self.ga['text'] = 'On'
            self.ga_config['state'] = 'normal'
        elif self.ga['text'] == 'On':
            self.ga['text'] = 'Off'
            self.ga_config['state'] = 'disabled'
            
            if 'num_generations' in self.configure:
                self.configure.pop('num_generations')
                self.configure.pop('num_parents_mating')
                self.configure.pop('sol_per_pop')
                self.configure.pop('pipelines')
                
                self.configure.pop('parent_selection_type')
                self.configure.pop('mutation_type')
                self.configure.pop('crossover_type')
                self.configure['split_half'] = False
                self.split['text'] = 'Off'
        
        if self.main_set_butt['state'] == 'normal':
            if 'num_generations' in self.configure or self.ga['text'] == 'Off':
                self.configure_now['state'] = 'normal'
                if self.data_dir != 'N/A' and self.output != 'N/A':
                    self.run_now['state'] = 'normal'
            else:
                self.configure_now['state'] = 'disabled'
                self.run_now['state'] = 'disabled'
            
    def genetic_config(self):
        self.configure['pipelines'] = int(self.pipelines.get())
        pipelines = self.to_num(self.pipelines.get())
        gen = 5
        sol = int(pipelines / gen)
        mate = int(sol / 2)
        
        self.genetic_frame_main = ttk.Frame(self.tab_main)
        self.genetic_frame_main.grid(row=0, column=0, sticky="nsew")
        self.genetic_frame_main.columnconfigure(0, weight=1)
        ttk.Label(self.genetic_frame_main, text='Genetic Algorithm Settings', font='Helvetica 12 bold').grid(row=0)
        
        self.genetic_frame = ttk.Frame(self.genetic_frame_main, relief='groove', padding='0.5i')
        self.genetic_frame.grid(row=1, column=0, sticky="nsew")
        self.genetic_frame.columnconfigure(0, weight=1)
        
        ttk.Label(self.genetic_frame, text='Number of generations:').grid(row=1)
        self.gen_num = ttk.Entry(self.genetic_frame)
        self.gen_num.insert(4, str(gen))
        self.gen_num.grid(row=1, column=1)
        
        ttk.Label(self.genetic_frame, text='Number of parents mating:').grid(row=2)
        self.p_num = ttk.Entry(self.genetic_frame)
        self.p_num.insert(4, str(mate))
        self.p_num.grid(row=2, column=1)
        
        ttk.Label(self.genetic_frame, text='Number of solutions per generation:').grid(row=3)
        self.sol_num = ttk.Entry(self.genetic_frame)
        self.sol_num.insert(4, str(sol))
        self.sol_num.grid(row=3, column=1)
        
        ttk.Label(self.genetic_frame, text='Parent selection type:').grid(row=4)
        selection_methods = ['sss', 'rws', 'sus', 'rank', 'random', 'tournament']
        self.parent_selection_var = tk.StringVar()
        self.parent_selection_var.set(selection_methods[0])
        self.parent_selection = tk.OptionMenu(self.genetic_frame, self.parent_selection_var, *selection_methods)
        self.parent_selection.grid(row=4, column=1)
        
        ttk.Label(self.genetic_frame, text='Crossover type:').grid(row=5)
        crossover_methods = ['single_point', 'two_points', 'uniform', 'scattered']
        self.cross_selection_var = tk.StringVar()
        self.cross_selection_var.set(crossover_methods[0])
        self.cross_selection = tk.OptionMenu(self.genetic_frame, self.cross_selection_var, *crossover_methods)
        self.cross_selection.grid(row=5, column=1)
        
        ttk.Label(self.genetic_frame, text='Mutation type:').grid(row=6)
        mutation_methods = ['random', 'swap', 'inversion', 'scramble', 'adaptive']
        self.mute_selection_var = tk.StringVar()
        self.mute_selection_var.set(mutation_methods[0])
        self.mute_selection = tk.OptionMenu(self.genetic_frame, self.mute_selection_var, *mutation_methods)
        self.mute_selection.grid(row=6, column=1)
        
        self.set_genetic = tk.Button(self.genetic_frame_main, text='Set', command=self.set_gen)
        self.set_genetic.grid(row=2)
        
        self.genetic_frame_main.tkraise()
        
    def set_gen(self):
        self.configure['num_generations'] = int(self.gen_num.get())
        self.configure['num_parents_mating'] = int(self.p_num.get())
        self.configure['sol_per_pop'] = int(self.sol_num.get())
        
        self.configure['parent_selection_type'] = self.parent_selection_var.get()
        self.configure['mutation_type'] = self.mute_selection_var.get()
        self.configure['crossover_type'] = self.cross_selection_var.get()
        
        self.configure['split_half'] = True
        self.split['text'] = 'On'
        
        if self.main_set_butt['state'] == 'normal':
            if 'num_generations' in self.configure or self.ga['text'] == 'Off':
                self.configure_now['state'] = 'normal'
                if self.data_dir != 'N/A' and self.output != 'N/A':
                    self.run_now['state'] = 'normal'
            else:
                self.configure_now['state'] = 'disabled'
                self.run_now['state'] = 'disabled'
        
        self.aesthetic_frame_main.tkraise()
        
    
    def define_seeds(self):
        self.define_frame = ttk.Frame(self.tab_main)
        self.define_frame.grid(row=0, column=0, sticky="nsew")
        self.define_frame.columnconfigure(0, weight=1)
        ttk.Label(self.define_frame, text='Seed Configuration', font='Helvetica 12 bold').grid(row=0)
        
        try:
            networks = int(self.networks.get())
        except:
            raise ValueError('Invalid number of networks. Please enter an integer - defaulting to 1 network')
            networks = 1
            
        modes = [self.atlas.get(), self.coords.get()]
        
        counter = 1
        
        self.out_dic['level1']['~construct~Finfo_rest_type']['gene'] = []
        self.out_dic['level1']['~construct~Finfo_rest_coords'] = {'gene': []}
        self.out_dic['level1']['~construct~Finfo_rest_seedinfo'] = {'gene': []}
        
        if self.atlas.get():
            self.out_dic['level1']['~construct~Finfo_rest_type']['gene'].append(self.atlas.get())
        if self.data.get():
            self.out_dic['level1']['~construct~Finfo_rest_type']['gene'].append(self.data.get())
            modes[0] = 1
        if self.coords.get():
            self.out_dic['level1']['~construct~Finfo_rest_type']['gene'].append(self.coords.get())
            
        if not sum(modes):
            raise ValueError('At least one seed definition type must be selected')
        
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
                        
                        vars(self)['atlas_add'+str(val)] = tk.Button(vars(self)['define_network'+str(val)], text='Add Brain Region', command=partial(self.add_option, i, val, 'called_'+str(i)+'_'+str(val), True, networks, modes))
                        vars(self)['atlas_add'+str(val)].grid(row=counter+val, column=1)
                        #counter += 1
                        
                        vars(self)['atlas_add'+str(val)] = tk.Button(vars(self)['define_network'+str(val)], text='Add Saved Mask', command=partial(self.add_option, i, val, 'called_'+str(i)+'_'+str(val), False, networks, modes))
                        vars(self)['atlas_add'+str(val)].grid(row=counter+val, column=0)
                        counter += 1
                    if i == 1:
                        vars(self)['atlas_dropdown'+str(val)] = []
                        fsl = os.getenv('FSLDIR')
                        for atlas in glob.glob(fsl + '/data/atlases/HarvardOxford*.xml'):
                            tree = ET.parse(atlas)
                            root = tree.getroot()
                        
                            for label in root.iter('label'):
                                vars(self)['atlas_dropdown'+str(val)].append(label.text)
                        
                        ttk.Label(vars(self)['define_network'+str(val)], text='MNI Voxel Coordinates (x, y, z):').grid(row=counter+val)
                        vars(self)['coords_'+str(val)] = tk.Entry(vars(self)['define_network'+str(val)])
                        vars(self)['coords_'+str(val)].insert(4, '')
                        vars(self)['coords_'+str(val)].grid(row=counter+val, column=1)
                        counter += 1
                        
                        vars(self)['called_'+str(i)+'_'+str(val)] = 0
                        ttk.Label(vars(self)['define_network'+str(val)], text='Multiverse Coordinate ROIs:').grid(row=counter+val)
                        vars(self)['called_'+str(i)+'_'+str(val)+'label'] = ttk.Label(vars(self)['define_network'+str(val)], text=str(vars(self)['called_'+str(i)+'_'+str(val)]))
                        vars(self)['called_'+str(i)+'_'+str(val)+'label'].grid(row=counter+val, column=1)
                        counter += 1
                        
                        vars(self)['coords_add'+str(val)] = tk.Button(vars(self)['define_network'+str(val)], text='Add Coordinates', command=partial(self.add_option, i, val, 'called_'+str(i)+'_'+str(val), True, networks, modes))
                        vars(self)['coords_add'+str(val)].grid(row=counter+val, column=1)
                        counter += 1
        
        self.main_set_butt = tk.Button(self.define_frame, text='Set', command=partial(self.ret, 'main'))
        self.main_set_butt['state'] = 'disabled'
        self.main_set_butt.grid(row=counter+val, column=0)
        
        self.define_frame.tkraise()
        
    def add_option(self, i, val, called_name, atlas, networks, modes):
        self.configure['networks'] = networks
        if i == 0:
            options = len(self.out_dic['level1']['~construct~Finfo_rest_seedinfo']) - 1
            brain_region = vars(self)['drop_select'+str(val)].get()
            if brain_region == '' or brain_region == ' ':
                raise ValueError('A brain region must be selected.')
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
                seedinfo = [tk.filedialog.askopenfilename()]
                dir = os.path.dirname(__file__)
                for l, seed in enumerate(seedinfo):
                    file_name = re.search('.*/(.*)', seed).group(1)
                    destination = os.path.join(dir, '.', 'seed_masks')
                    if not os.path.isdir(destination):
                        os.mkdir(destination)
                    
                    copy2(seed, os.path.join(destination, file_name))
                    seedinfo[l] = os.path.join('.', 'seed_masks', file_name)
                    
                seedinfo = seedinfo[0]
                
                
            if vars(self)[called_name] < options:
                for j in range(vars(self)[called_name], options):
                    if len(self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][j]) <= val:
                        self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][j].append(seedinfo)
                    else:
                        self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][j][val] = seedinfo
            else:
                if options > 0:
                    self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options] = self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options-1].copy()
                    self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options][val] = seedinfo
                else:
                    self.out_dic['level1']['~construct~Finfo_rest_seedinfo'][options] = [seedinfo] * networks

            self.out_dic['level1']['~construct~Finfo_rest_seedinfo']['gene'] = {'low': -0.5, 'high': options+0.5}
            vars(self)['drop_select'+str(val)].set(" ")
            vars(self)['atlas_max'+str(val)].delete(0, 'end')
            vars(self)['atlas_max'+str(val)].insert(4, '95')
            vars(self)['atlas_min'+str(val)].delete(0, 'end')
            vars(self)['atlas_min'+str(val)].insert(4, '0')
            
        elif i == 1:
            options = len(self.out_dic['level1']['~construct~Finfo_rest_coords']) - 1
            coords = vars(self)['coords_'+str(val)].get()
            coords = coords.replace(' ', '').split(',')
            if len(coords) != 3:
                raise ValueError('3 values must be entered as a comma separated list corresponding to the x, y, and z axis')
                
            for i, k in enumerate(coords):
                coords[i] = self.to_num(k)
                
            coords = tuple(coords)
            
            if vars(self)[called_name] < options:
                for j in range(vars(self)[called_name], options):
                    if len(self.out_dic['level1']['~construct~Finfo_rest_coords'][j]) <= val:
                        self.out_dic['level1']['~construct~Finfo_rest_coords'][j].append(coords)
                    else:
                        self.out_dic['level1']['~construct~Finfo_rest_coords'][j][val] = coords
            else:
                if options > 0:
                    self.out_dic['level1']['~construct~Finfo_rest_coords'][options] = self.out_dic['level1']['~construct~Finfo_rest_coords'][options-1].copy()
                    self.out_dic['level1']['~construct~Finfo_rest_coords'][options][val] = coords
                else:
                    self.out_dic['level1']['~construct~Finfo_rest_coords'][options] = [coords] * networks
            
            self.out_dic['level1']['~construct~Finfo_rest_coords']['gene'] = list(range(len(self.out_dic['level1']['~construct~Finfo_rest_coords'])-1))
            
            vars(self)['coords_'+str(val)].delete(0, 'end')
            vars(self)['coords_'+str(val)].insert(8, '')
            
        vars(self)[called_name] += 1
        vars(self)[called_name+'label']['text'] = vars(self)[called_name]
        
        check = True
        for val in range(networks):
            for i, p in enumerate(modes):
                if p:
                    if not vars(self)['called_'+str(i)+'_'+str(val)]:
                        check = False
        
        if check:
            self.main_set_butt['state'] = 'normal'
            if 'num_generations' in self.configure or self.ga['text'] == 'Off':
                self.configure_now['state'] = 'normal'
                if self.data_dir != 'N/A' and self.output != 'N/A':
                    self.run_now['state'] = 'normal'

    
    def get_format(self):
        dir = os.path.dirname(__file__)
        #NOTE: CONFIG FILE OUTSIDE OF FOLDER, PUT CONFIGURATION IN CODE
        path = os.path.join(dir, '../configuration/default.json')
        with open(path) as f:
            data = json.load(f)
        
        return data
    
    def strings(self, stage, alias, default, value_map, param_name_full, value, butt_text, counter):
        vars(self)['aesthetic_frame_'+stage+param_name_full] = ttk.Frame(vars(self)['tab_window_'+stage], relief="groove", padding="0.5i")
        vars(self)['aesthetic_frame_'+stage+param_name_full].grid(row=counter, column=0, sticky="nsew")
        vars(self)['aesthetic_frame_'+stage+param_name_full].columnconfigure(0, weight=1)
        
        ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text=alias).grid(row=1)
        
        if butt_text == "Off":
            text = str(value_map[default])
        else:
            text = str(value_map)
            
        ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text='Value Pool:').grid(row=2)
        vars(self)['valuepool_'+stage+param_name_full] = ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text=text)
        vars(self)['valuepool_'+stage+param_name_full].grid(row=2, column=1)

        vars(self)[param_name_full+stage+'configbutton'] = tk.Button(vars(self)['aesthetic_frame_'+stage+param_name_full], text='Configure', 
                                                                command=partial(self.mapped_config_frame, stage, param_name_full, alias, value, value_map))
        vars(self)[param_name_full+stage+'configbutton'].grid(row=3, column=1)
        
        vars(self)[param_name_full+stage+'onbutton'] = tk.Button(vars(self)['aesthetic_frame_'+stage+param_name_full], text=butt_text, command=partial(self.change, param_name_full+stage+'onbutton', 'valuepool_'+stage+param_name_full, default, value, value_map, stage, param_name_full, param_name_full+stage+'configbutton'))
        vars(self)[param_name_full+stage+'onbutton'].grid(row=1, column=1)
        
        
    
    def mapped_config_frame(self, stage, param_name_full, alias, value, value_map):
        vars(self)[param_name_full+stage+'frame'] = tk.Frame(vars(self)['tab_'+stage])
        vars(self)[param_name_full+stage+'frame'].grid(row=0, column=0, sticky="nsew")
        vars(self)[param_name_full+stage+'frame'].columnconfigure(0, weight=1)
        ttk.Label(vars(self)[param_name_full+stage+'frame'], text=alias + ' Configuration', font='Helvetica 12 bold').grid(row=0)
        
        for val in value:
            vars(self)[param_name_full+stage+str(val)+'_var'] = tk.IntVar()
            vars(self)[param_name_full+stage+str(val)] = tk.Checkbutton(vars(self)[param_name_full+stage+'frame'], text=str(value_map[val]), variable=vars(self)[param_name_full+stage+str(val)+'_var'])
            vars(self)[param_name_full+stage+str(val)].grid(row=val+1)
            
        vars(self)[param_name_full+stage+'frame_button'] = tk.Button(vars(self)[param_name_full+stage+'frame'], text='Exit', 
                                                                   command=partial(self.ret, stage))
        vars(self)[param_name_full+stage+'frame_button'].grid(row=val+2, column=0)
        
        vars(self)[param_name_full+stage+'frame_button'] = tk.Button(vars(self)[param_name_full+stage+'frame'], text='Set', 
                                                               command=partial(self.selection, value, value_map, param_name_full, stage))
        vars(self)[param_name_full+stage+'frame_button'].grid(row=val+2, column=1)
        
        vars(self)[param_name_full+stage+'frame'].tkraise()
        
            
    def selection(self, value, value_map, param_name_full, stage):
        lst = []
        updated = []
        for val in value:
            temp = vars(self)[param_name_full+stage+str(val)+'_var'].get()
            if temp:
                lst.append(val)
                updated.append(value_map[val])
        if not lst:
            raise ValueError('No values selected. Must select 1 or more')
            
        self.out_dic[stage][param_name_full]['gene'] = lst
        vars(self)['valuepool_'+stage+param_name_full]['text'] = str(updated)
        
        vars(self)['tab_window_'+stage].tkraise()
        
    
    def numerical(self, stage, alias, default, value_map, param_name_full, value, butt_text, counter):
        vars(self)['aesthetic_frame_'+stage+param_name_full] = ttk.Frame(vars(self)['tab_window_'+stage], relief="groove", padding="0.5i")
        vars(self)['aesthetic_frame_'+stage+param_name_full].grid(row=counter, column=0, sticky="nsew")
        vars(self)['aesthetic_frame_'+stage+param_name_full].columnconfigure(0, weight=1)
        
        ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text=alias).grid(row=1)
        
        ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text='Value Pool:').grid(row=2)
        
        if butt_text == "Off":
            text = str(value[default])
        else:
            text = str(value)
            
        vars(self)['valuepool_'+stage+param_name_full] = ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text=text)
        vars(self)['valuepool_'+stage+param_name_full].grid(row=2, column=1)

        ttk.Label(vars(self)['aesthetic_frame_'+stage+param_name_full], text='Input Type:').grid(row=3)
        vars(self)[param_name_full+stage+'rangebutton'] = tk.Button(vars(self)['aesthetic_frame_'+stage+param_name_full], text='Range', command=partial(self.change, param_name_full+stage+'rangebutton', 'valuepool_'+stage+param_name_full, default, value, value_map, stage, param_name_full, ''))
        vars(self)[param_name_full+stage+'rangebutton'].grid(row=3, column=1)
        
        vars(self)[param_name_full+stage+'configbutton'] = tk.Button(vars(self)['aesthetic_frame_'+stage+param_name_full], text='Configure', 
                                                                command=partial(self.config_frame, param_name_full+stage+'rangebutton', stage, param_name_full, alias))
        vars(self)[param_name_full+stage+'configbutton'].grid(row=4, column=1)
        
        vars(self)[param_name_full+stage+'onbutton'] = tk.Button(vars(self)['aesthetic_frame_'+stage+param_name_full], text=butt_text, command=partial(self.change, param_name_full+stage+'onbutton', 'valuepool_'+stage+param_name_full, default, value, value_map, stage, param_name_full, param_name_full+stage+'configbutton'))
        vars(self)[param_name_full+stage+'onbutton'].grid(row=1, column=1)

    def change(self, name, def_name, default, value, value_map, stage, param_name_full, configure):
        if vars(self)[name]['text'] == 'On':
            vars(self)[name].config(text='Off')
            vars(self)[configure]['state'] = 'disabled'
            if value_map:
                vars(self)[def_name].config(text=str(value_map[default]))
                self.out_dic[stage][param_name_full]['gene'] = [default]
            else:
                vars(self)[def_name].config(text=str(value[default]))
                self.out_dic[stage][param_name_full]['gene'] = [value[default]]
        elif vars(self)[name]['text'] == 'Off':
            vars(self)[name].config(text='On')
            vars(self)[configure]['state'] = 'normal'
            if value_map:
                vars(self)[def_name].config(text=str(value_map))
                self.out_dic[stage][param_name_full]['gene'] = value
            else:
                vars(self)[def_name].config(text=str(value))
                self.out_dic[stage][param_name_full]['gene'] = value
        elif vars(self)[name]['text'] == 'Range':
            vars(self)[name].config(text='Manual')
        elif vars(self)[name]['text'] == 'Manual':
            vars(self)[name].config(text='Range')
                
    def config_frame(self, range_state, stage, param_name_full, alias):
        vars(self)[param_name_full+stage+'frame'] = tk.Frame(vars(self)['tab_'+stage])
        vars(self)[param_name_full+stage+'frame'].grid(row=0, column=0, sticky="nsew")
        vars(self)[param_name_full+stage+'frame'].columnconfigure(0, weight=1)
        ttk.Label(vars(self)[param_name_full+stage+'frame'], text=alias + ' Configuration', font='Helvetica 12 bold').grid(row=0)
        if vars(self)[range_state]['text'] == 'Range':
            ttk.Label(vars(self)[param_name_full+stage+'frame'], text='Low: ').grid(row=1)
            vars(self)[param_name_full+stage+'low'] = tk.Entry(vars(self)[param_name_full+stage+'frame'])
            vars(self)[param_name_full+stage+'low'].grid(row=1, column=1)
            ttk.Label(vars(self)[param_name_full+stage+'frame'], text='High: ').grid(row=2)
            vars(self)[param_name_full+stage+'high'] = tk.Entry(vars(self)[param_name_full+stage+'frame'])
            vars(self)[param_name_full+stage+'high'].grid(row=2, column=1)
            ttk.Label(vars(self)[param_name_full+stage+'frame'], text='Step: ').grid(row=3)
            vars(self)[param_name_full+stage+'step'] = tk.Entry(vars(self)[param_name_full+stage+'frame'])
            vars(self)[param_name_full+stage+'step'].grid(row=3, column=1)
            
            vars(self)[param_name_full+stage+'frame_button'] = tk.Button(vars(self)[param_name_full+stage+'frame'], text='Exit', 
                                                                   command=partial(self.ret, stage))
            vars(self)[param_name_full+stage+'frame_button'].grid(row=4, column=0)
            
            vars(self)[param_name_full+stage+'frame_button'] = tk.Button(vars(self)[param_name_full+stage+'frame'], text='Set', 
                                                                   command=partial(self.set_values, range_state, param_name_full, stage))
            vars(self)[param_name_full+stage+'frame_button'].grid(row=4, column=1)
        elif vars(self)[range_state]['text'] == 'Manual':
            ttk.Label(vars(self)[param_name_full+stage+'frame'], text='Values (comma separated list): ').grid(row=1)
            vars(self)[param_name_full+stage+'manual'] = tk.Entry(vars(self)[param_name_full+stage+'frame'])
            vars(self)[param_name_full+stage+'manual'].grid(row=1, column=1)
            
            vars(self)[param_name_full+stage+'frame_button'] = tk.Button(vars(self)[param_name_full+stage+'frame'], text='Exit', 
                                                                   command=partial(self.ret, stage))
            vars(self)[param_name_full+stage+'frame_button'].grid(row=2, column=0)
            
            vars(self)[param_name_full+stage+'frame_button'] = tk.Button(vars(self)[param_name_full+stage+'frame'], text='Set', 
                                                                   command=partial(self.set_values, range_state, param_name_full, stage))
            vars(self)[param_name_full+stage+'frame_button'].grid(row=2, column=1)
            
        vars(self)[param_name_full+stage+'frame'].tkraise()
            
    def set_values(self, range_state, param_name_full, stage):
        if vars(self)[range_state]['text'] == 'Range':
            low = self.to_num(vars(self)[param_name_full+stage+'low'].get())
            high = self.to_num(vars(self)[param_name_full+stage+'high'].get())
            step = self.to_num(vars(self)[param_name_full+stage+'step'].get())
            if step:
                self.out_dic[stage][param_name_full]['gene'] = {'low': low, 'high': high, 'step': step}
                updated = {'low': low, 'high': high, 'step': step}
            else:
                self.out_dic[stage][param_name_full]['gene'] = {'low': low, 'high': high}
                updated = {'low': low, 'high': high}
            
            vars(self)['valuepool_'+stage+param_name_full]['text'] = str(step)
        elif vars(self)[range_state]['text'] == 'Manual':
            updated = vars(self)[param_name_full+stage+'manual'].get()
            updated = updated.replace(' ', '').split(',')
            for i, val in enumerate(updated):
                updated[i] = self.to_num(val)
            
            self.out_dic[stage][param_name_full]['gene'] = updated
            
        vars(self)['valuepool_'+stage+param_name_full]['text'] = str(updated)
        vars(self)['tab_window_'+stage].tkraise()
    
    def ret(self, stage):
        if stage == 'main':
            self.aesthetic_frame_main.tkraise()
        else:
            vars(self)['tab_window_'+stage].tkraise()
            
    def to_num(self, str_):
        try:
            return int(str_)
        except:
            try:
                return float(str_)
            except:
                raise ValueError('Invalid input type, must be int or float')
                