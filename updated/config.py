#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:45:59 2021

@author: grahamseasons
"""
# =============================================================================
import subprocess
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

create_docker_command = '-v {code}:/root/multiverse/code -v {data}:/data -v {work_dir}:/scratch -w scratch'
create_docker_command = create_docker_command.format(code=code, data=data_dir, work_dir=exp_dir)
#THIS SHOULD SET UP DOCKER, AND ANALYSIS CODE, DOCKER SHOULD THEN AUTOMATICALLY RUN ANALYSIS CODE
subprocess.run([])