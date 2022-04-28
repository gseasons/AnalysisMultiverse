#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:12:07 2022

@author: grahamseasons
"""
import sys
from functions import organize
import pickle
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import nibabel as nib
import numpy as np

processed = '/scratch/processed'
    
if len(sys.argv) > 2:
    task = sys.argv[1]
    sq = sys.argv[2]
    
#out_frame = processed + '/{task}.pkl'.format(task=task)
out_frame = '/Volumes/NewVolume/_i_0/rest.pkl'
task= 'rest'
if True:#sq.count('batch.sh') == 1:
    paths = organize(task, out_frame)
    
    with open(paths, 'rb') as f:
        paths = pickle.load(f)
    #HARDCODED
    with open('/Volumes/NewVolume/_i_0/generation_0.pkl', 'rb') as f:
        pipelines = pickle.load(f)
        
    mds = MDS()
    embedded = mds.fit_transform(pipelines.transpose())
    fig = plt.figure(1)
    ax = plt.axes()
    plt.scatter(embedded[:,0], embedded[:,1])
    plt.scatter(embedded[:8,0], embedded[:8,1])
    
    segments = []
    similarity = []
    
    for i, x in enumerate(embedded[:8]):
        for j, y in enumerate(embedded[:8]):
            if i >= j:
                continue
            
            paths1 = paths['pipeline'][i]['network']
            out_stats = []
            #
            for k in range(len(paths['pipeline'][i]['network'])):
                for l, m in enumerate(list(paths['pipeline'][i]['network'][k]['contrast'].values())):
                    m2 = list(paths['pipeline'][j]['network'][k]['contrast'].values())[l]
                    out = np.corrcoef(np.array(nib.load(m).dataobj).flatten(), np.array(nib.load(m2).dataobj).flatten())[0,1]
                    out_stats.append(out)
            
            out = [abs(stat) for stat in out_stats if not np.isnan(stat)]
            out = np.mean(out_stats)
            if np.isnan(out):
                out = 0
                continue
            else:
                pass
                #out = abs(out)
                
            similarity.append(out)
            
            segments.append([list(x), list(y)])
        
    similarity = np.array(similarity)
    values = np.abs(similarity)
    lc = LineCollection(segments, linewidths=np.full(len(segments), 1), zorder=0, cmap=plt.cm.Purples, norm=plt.Normalize(0, values.max()))
    lc.set_array(similarity)
    ax.add_collection(lc)
    
    plt.show()
    A=3
            
        
    
    
    #TODO: THE THIRD # POINT BELOW THIS ONE
    
    #INSERT DATA PROCESSING HERE - TO BE WRITTEN ONCE WE HAVE TEST DATA
    
    #MULTIDIMENSIONAL SCALING TO REDUCE DIMENSIONS FROM ~60 -> 2 ( or t-SNE )
    #FEED IN GA MATRIX -> REPLACE DEPENDENT VARIABLES (i.e. ROI variables if it's REHO) WITH DEFAULT VALUES (FOR REST)
    
    
    #CLUSTERS PIPELINES IN 2D SO CAN SEE SIMILARITY IN TERMS OF CHOSEN PARAMETERS
    #FOR DATA PRESENTATION COULD ADD ANOTHER DIMENSION OF SIMILARITY SCORE WITH PIPELINES IN NEIGHBOURHOOD (SIMILARITY WEIGHTED BY EUCLIDEAN DISTANCE TO EACH POINT)
    #POOR SIMILARITY IN NEIGHBOURHOOD, CATEGORIZE DIFFERENCES BETWEEN PIPELINES TO COMPARE
    #FIND ACTUAL BEST SIMILARITY BETWEEN PIPELINES (CORRELATION MATRIX TYPE THING?)
    #SHOW SIMILARITY ON PLOT AS LINE THICKNESS/OPACITY https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#sphx-glr-auto-examples-manifold-plot-mds-py
    #CHECK TO SEE WHICH PARAMETERS DIFFER BETWEEN POORLY PERFORMING PIPELINES IN SAME GROUP -> REORGANIZE/REPLOT ONLY CONSIDERING THOSE
    #   SAME IDEA WITH ONES THAT DO MATCH WELL
    
    #OPTIMIZATION ALGORITHM FOR ABOVE?
    
    #SCATTER PLOT FOR EACH BRAIN REGION WE GET SIGNIFICANT ACTIVATION -> SHOW WHICH PIPELINES AGREE (others in faded colour)
    
    #REORGANIZE ORIGINAL MDS BY SIMILARITY OF OUTPUTS
    #SHOW AVERAGE MAPS (or CORRELATION MATRICES OF BRAIN REGIONS) OR AN INTERSECTION OF ACTIVATION MAPS
    
    #IF NONE OF THEM ARE SPATIALLY SIMILAR, OR ACTIVATION MAPS AREN'T SIMILAR: TBD
    #STAINED GLASS AVERAGE (ALL DIFERENT COLOURS, TRANSPARENCIES OVERLAYED)?
    #MAP OF BRAIN DISPLAYING AVERAGE ACTIVATION ACROSS PIPELINES, AND MAP WHERE EACH VOXEL REPRESENTS THE VARIANCE ACROSS PIPELINES