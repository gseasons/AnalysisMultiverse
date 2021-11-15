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
#SAVE CONFIG FILE WITH DEFAULTS AND ALTERED VALUES TO LOCATION
subprocess.run([])

# =============================================================================
# neurodocker generate docker --base debian:stretch --pkg-manager apt --fsl version=6.0.0 --ants version=2.3.1 --afni version=latest --user multiverse --miniconda miniconda_version="4.3.31" create_env=multiverse conda_install="python=3.8 pytest traits pandas matplotlib scikit-learn scikit-image seaborn numpy" pip_install="nipype nilearn nibabel niworkflows pygad kivy" --user root --install nano --volume /data /root/multiverse/code /scratch -w /scratch > multiverse.Dockerfile
# docker build --tag multiverse --file multiverse.Dockerfile .
# =============================================================================
docker run -v /Users/grahamseasons/fMRI/analysis_multiverse/updated:/root/multiverse/code -v /Volumes/NewVolume/super_agers:/data -v /Volumes/NewVolume/docker_test:/scratch -w /scratch multiverse