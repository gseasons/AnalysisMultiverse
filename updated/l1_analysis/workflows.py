#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:50:42 2021

@author: grahamseasons
"""

#have so that multiple values inputted as tuples
#resting will be a dictionary with everything seed related inside
#multiple seeds, separate statistical maps

def sessioninfo(mask, task, TR, event_file, unsmoothed, smoothed, outliers, brain, brainmask, segmentations, invwarp, network):
    from nipype import Node
    from versatile import SpecifyModelVersatile
    from nipype.interfaces.fsl import ImageMeants, ExtractROI, ImageMaths, ApplyWarp, Threshold, BinaryMaths, WarpPointsFromStd
    from nipype.interfaces.fsl.maths import MathsCommand
    from functions import parse_xml
    import re, os
    import numpy as np
    from updated.l1_analysis.functions import data_driven, warp, invert
    
    model = Node(SpecifyModelVersatile(input_units='secs', parameter_source='FSL'), name='model')
    model.inputs.outlier_files = outliers
    model.inputs.time_repetition = TR
    model.inputs.high_pass_filter_cutoff = vars().get('HP', 128)
    model.inputs.functional_runs = smoothed
    
    if event_file:
        model.inputs.bids_event_file = event_file
        session_info = model.run().outputs.session_info
    elif 'rest' in task and vars().get('rest'):
        #NOTE: COORDINATES WILL BE IN STANDARD SPACE, AS WILL ATLAS LOCATIONS -> TRANSFER TO STANDARD TO GET TIME COURSES
        if vars()['rest']['type'] == 'ROI':
            x, y, z = vars()['rest']['coords'][network]
            
            if vars()['warppostfeat']:
                roipoints = os.getcwd() + '/roipoints.txt'
                np.savetxt(roipoints, np.array(x,y,z))
                warppoints = WarpPointsFromStd(in_coords=roipoints, img_file=brain, std_file=mask, warp_file=invwarp).run().outputs.out_file
                points = np.loadtxt(warppoints)
                x = points[0]
                y = points[1]
                z = points[2]
            
            radius = vars()['rest']['radius']
            createseed = Node(MathsCommand(in_file=mask), name='createseed')
            createseed.inputs.args = '-mul 0 -add 1 -roi {x} 1 {y} 1 {z} 1 0 1'.format(x=x, y=y, z=z)
            seed = createseed.run().outputs.out_file
            
            makesphere = Node(MathsCommand(in_file=seed), name='makesphere')
            makesphere.inputs.args = '-kernel sphere {radius} -fmean'.format(radius=radius)
            
            thrseed = Node(ImageMaths(op_string='-bin'), name='thrseed')
            thrseed.inputs.in_file = makesphere.run().outputs.out_file
            thrfile = thrseed.run().outputs.out_file
            suffix = '_x{x}_y{y}_z{z}_r{r}'.format(x=x, y=y, z=z, r=radius)
        elif vars()['rest']['type'] == 'atlas' or vars()['rest']['type'] == 'data':
            #nipype.algorithms.misc PickAtlas(although i've kind of implemented this, might be nice to add multiple ROIs -> get user input for options and index them in GA, pass list of list)
            #ADD WHEN SETTING UP USER INPUT GUI -> atlas will be link to atlas to be used by PickAtlas
            seed, thr = vars()['rest']['seedinfo'][network]
            atlas = vars()['rest']['atlas']
            file, index, name = parse_xml(atlas, seed, mask)
            getseed = Node(ExtractROI(in_file=file, t_min=int(index), t_size=1), name='getseed')
            roi = getseed.run().outputs.roi_file
            
            if vars()['warppostfeat']:
                roi = ApplyWarp(in_file=roi, ref_file=brain, field_file=invwarp).run().outputs.out_file
            
            thrseed = Node(ImageMaths(op_string='-thr {thr} -bin'.format(thr=thr)), name='thrseed')
            thrseed.inputs.in_file = roi
            thrfile = thrseed.run().outputs.out_file
            suffix = ''.join([word[0] for word in name.split()])
            
            if vars()['rest']['type'] == 'data':
                k = vars()['rest']['k']
                kcc = vars()['rest']['kcc']
                lp = vars()['rest'].get('lp', 0.01)
                hp = vars()['rest'].get('hp', 0.1)
                
                reho = data_driven(brainmask, unsmoothed, k, kcc, TR, lp, hp)
                mul = Node(BinaryMaths(in_file=reho, operand_file=thrfile, operation='mul'), name='mul')
                thrfile = mul.run().outputs.out_file
                suffix += '_reho'
        else:
            print('ERROR')
            
    else:
        print('ERROR')
                
    return session_info

