#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 10:33:33 2021

@author: grahamseasons
"""

import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from os.path import abspath


def plot_slice(fname):

    # Load the image
    img = nb.load(fname)
    data = img.get_data()

    # Cut in the middle of the brain
    cut = int(data.shape[-1]/2) + 10

    # Plot the data
    plt.imshow(np.rot90(data[..., cut]), cmap="gray")
    plt.gca().set_axis_off()



#INTERFACES
from nipype.interfaces import fsl
from nipype.interfaces.fsl import BET
from nipype.interfaces.fsl import IsotropicSmooth
from nipype import Node, Function, Workflow

#NODES

bet = Node(BET(frac=0.3), name='bet_node')
bet.inputs.in_file = ""
bet.inputs.out_file = ""

smooth_node = Node(IsotropicSmooth(), name="smoothing")
smooth_node.inputs.in_file = ""
smooth_node.inputs.fwhm = 4
smooth_node.inputs.out_file = ""

#WORKFLOW 
# passing interfaces between each other:
# don't name intermediate files, done behind scenes, pass example.result.outputs.TYPE_file (example = node.run())
# into in_file of subsequent node.  Passing objects!
#
# Workflow
# file names as absolute paths to node important
# allows you to connect nodes and only specify the generic input/output categories
# need to specify a base directory which is where output images of workflow will be saved (under subfolders from the named nodes)

##### DATA GRABBING #####

# mandate bids format
from bids.layout import BIDSLayout
from nipype.interfaces.io import BIDSDataGrabber
#WILL WANT TO WRAP IN A NODE
data_struct = BIDSLayout("data/TBD") #Grab input data folder from user
subj = data_struct.get_subjects()
mod = data_struct.get_modalities()
im_type = data_struct.get_types(modality='func') #change modality depending on info for anat, dwi, func
tasks = data_struct.get_tasks()

# CAN ITERATE OVER ALL SUBJECTS & TASKS
bg_all = Node(BIDSDataGrabber(), name='bids-grabber')
bg_all.inputs.base_dir = '/data/ds000114'
bg_all.inputs.output_query = {'bolds': dict(type='bold', taskname=tasks[0])} #can also grab T1W images for structural
bg_all.iterables = ('subject', data_struct.get_subjects()[:])

bg_all.bolds #'Contains a FILE LIST OF ALL BOLD FILES'
met_dat = data_struct.get_metadata(bg_all.bolds[0])
TE = met_dat['EchoTime']
flip = met_dat['FlipAngle']
TR = met_dat['RepetitionTime']
ST = met_dat['SliceTiming']

#CALL PREPROCESS 
#SHOULD BE A WAY TO GET SESSION INFO AS WELL


##### PREPROCESS #####

# input and define parameters which define workflow #
fwhm = 4
iso_resample = 4
#TR = 5
HP = 50 
sig_loss_thresh = 20 #SIGNAL LOSS THRESHOLD 
bright_thresh = 2 #greater than noise, less than contrast of edges -> ask erin how she chooses/calculates


# end inputs #


from nipype.interfaces.fsl import (BET, ExtractROI, FAST, FLIRT, ImageMaths,
                                   MCFLIRT, SliceTimer, Threshold)
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.algorithms.rapidart import ArtifactDetect
from nipype import Workflow, Node, MapNode, JoinNode
from nipype.interfaces.fsl import IsotropicSmooth, SUSAN

exp_dir = '/output' #TO BE GRABBED FROM USER, base directory where all analysis conducted
out_dir = 'datasink' #sub folder inside of exp_dir, output of datasink
work_dir = 'workingdir' #where all intermediate files will be stored

#skip dummy scans, change t_min
extract = Node(ExtractROI(t_min=4, t_size=-1, output_type='NIFTI'), name='extract')
#motion correction
mc = Node(MCFLIRT(mean_vol=True, save_plots=True, output_type='NIFTI'), name='mc')
#slice timing correction
st = Node(SliceTimer(index_dir=False, interleaved=True, output_type='NIFTI', time_repetition=TR), name='st')
#smoothing (i think isotropic smooth is an option as well)
susan = Node(SUSAN(brightness_threshold=bright_thresh, fwhm=fwhm), name='susan')
#Artifact detection - not fsl native but seems good, introduces more parameters though
art = Node(ArtifactDetect(norm_threshold=2,
                          zintensity_threshold=3,
                          mask_type='spm_global',
                          parameter_source='FSL',
                          use_differences=[True, False],
                          plot_type='svg'),
           name="art")



#coregistration workflow


















