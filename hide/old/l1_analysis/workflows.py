#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:50:42 2021

@author: grahamseasons
"""

def info(mask, task, TR, event_file, unsmoothed, smoothed, brain, brainmask, outliers, segmentations, invwarp, network):
    from nipype import Node
    from versatile import SpecifyModelVersatile
    from nipype.interfaces.fsl import ImageMeants, ExtractROI, ImageMaths, BinaryMaths, WarpPoints
    from nipype.interfaces.fsl.maths import MathsCommand
    from l1_analysis.functions import parse_xml
    import re, os
    import numpy as np
    from l1_analysis.functions import data_driven, warp
    
    model = Node(SpecifyModelVersatile(input_units='secs', parameter_source='FSL'), name='model')
    model.inputs.time_repetition = TR
    model.inputs.high_pass_filter_cutoff = vars().get('HP', 128)
    model.inputs.functional_runs = smoothed
    model.inputs.outlier_files = outliers
    
    if event_file:
        model.inputs.bids_event_file = event_file
        session_info = model.run().outputs.session_info
        session_info[0]['scans'] = smoothed
        return session_info
    elif 'rest' in task and vars().get('rest'):
        if vars()['rest']['type'] == 'ROI':
            x, y, z = vars()['rest']['coords'][network]
            
            if vars()['warppostfeat']:
                roipoints = os.getcwd() + '/roipoints.txt'
                np.savetxt(roipoints, np.array([x,y,z]))
                if vars()['concatenate']:
                    img_file = brainmask
                    mask_ = brainmask
                else:
                    img_file = brain
                    mask_ = brain
                warppoints = WarpPoints(in_coords=roipoints, dest_file=img_file, src_file=mask, warp_file=invwarp).run().outputs.out_file
                points = np.loadtxt(warppoints).astype(int)
                x = points[0]
                y = points[1]
                z = points[2]
            else:
                mask_ = mask
            
            radius = vars()['rest']['radius']
            createseed = Node(MathsCommand(in_file=mask_, output_datatype='float'), name='createseed')
            createseed.inputs.args = '-mul 0 -add 1 -roi {x} 1 {y} 1 {z} 1 0 1'.format(x=x, y=y, z=z)
            seed = createseed.run().outputs.out_file
            
            makesphere = Node(MathsCommand(in_file=seed), name='makesphere')
            makesphere.inputs.args = '-kernel sphere {radius} -fmean'.format(radius=radius)
            
            thrseed = ImageMaths(op_string='-bin')
            thrseed.inputs.in_file = makesphere.run().outputs.out_file
            thrfile = thrseed.run().outputs.out_file
            suffix = '_x{x}_y{y}_z{z}_r{r}'.format(x=x, y=y, z=z, r=radius)
        elif vars()['rest']['type'] == 'atlas' or vars()['rest']['type'] == 'data':
            #nipype.algorithms.misc PickAtlas(although i've kind of implemented this, might be nice to add multiple ROIs -> get user input for options and index them in GA, pass list of list)
            #ADD WHEN SETTING UP USER INPUT GUI -> atlas will be link to atlas to be used by PickAtlas
            if isinstance(vars()['rest']['seedinfo'][network], (tuple, list)):
                seed, thr = vars()['rest']['seedinfo'][network]
                atlas = vars()['rest']['atlas']
                file, index, name = parse_xml(atlas, seed, mask)
                suffix = ''.join([word[0] for word in name.split()])
                getseed = Node(ExtractROI(in_file=file, t_min=int(index), t_size=1), name='getseed')
                roi = getseed.run().outputs.roi_file
            else:
                roi = vars()['rest']['seedinfo'][network]
                suffix = re.search('.*([A-Za-z0-9_-]+).nii.*').group(1)
                thr = 0
            
            if vars()['warppostfeat']:
                if vars()['concatenate']:
                    ref_file = brainmask
                else:
                    ref_file = brain
                roi = warp(roi, ref_file, invwarp)
            else:
                ref_file = brainmask
            
            thrseed = ImageMaths(op_string='-thr {thr} -bin'.format(thr=thr))
            thrseed.inputs.in_file = roi
            thrfile = thrseed.run().outputs.out_file
            
            
            if vars()['rest']['type'] == 'data':
                k = vars()['rest']['k']
                kcc = vars()['rest']['kcc']
                lp = vars()['rest'].get('lp', 0.01)
                hp = vars()['rest'].get('hp', 0.1)
                
                reho = data_driven(ref_file, unsmoothed, k, kcc, TR, lp, hp)
                mul = BinaryMaths(in_file=reho, operand_file=thrfile, operation='mul')
                thrfile = mul.run().outputs.out_file
                suffix += '_reho'
        else:
            raise NotImplementedError("Invalid seeding method of {method} used, which is not implemented. Please use 'atlas', 'ROI', or 'data'.".format(method=vars()['rest']['type']))
            
        ev_name = re.search('task-([a-zA-Z]+)_', smoothed).group(1) + suffix
                    
        mean_ts = Node(ImageMeants(in_file=smoothed, out_file=ev_name), name='mean_ts')
                    
        mean_ts.inputs.mask = thrfile
        time_series = [mean_ts.run().outputs.out_file]
                
        model.inputs.event_files = time_series
        session_info = model.run().outputs.session_info
        
        session_info[0]['scans'] = smoothed
        
        return session_info, thrfile
    else:
        raise ValueError("Unhandled task of {task}. If resting state analysis ensure 'rest' is in the task name, otherwise ensure there is a valid event_file".format(task=task))
    
