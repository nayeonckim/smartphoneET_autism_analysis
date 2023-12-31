#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute gaze time to bodyparts obtained from densepose.
In the current analysis, this code was used to identify the 'whole body' areas.
Areas in each movie frame that were not included in the 'whole body' (i.e. 'non-body')
were annotated as 'background' in the manuscript 

adapted and modified from https://github.com/adolphslab/adolphslab_eyetracking

"""
import sys
import os
import numpy as np
import pickle
from tqdm import tqdm

# The Adolphs Lab data analysis scripts for eye tracking are provided under alabeye package  
from alabeye.etdata import ETdata, makedir
from alabeye.io.loadstimfeatures import densepose_results2bodypart_layers, densepose_facehand_layers, radmask

# if matlab version of output is needed
from scipy.io import savemat 

# Main directory for experiment data
root_dir = 'path-to-dir'

# The main directory for experimental data contains several subdirectories such as:
# - ETdata: eye tracking data
# - ETdata/down2frame_data: gaze data downsampled to the frame rate of video stimulus
# - BehavioralData: demographic and psychological assessment info about subjects
# - StimVids: media files for experimental stimuli [not shared here because of copyright restrictions] 
# - FeatureZoo: various features extracted from video stimuli 

stim_dir = 'path-to-stim'
features_dir = 'path-to-features' # feature annotation saved in this directory
prepdata_dir = os.path.join(root_dir,'tobii','down2frame_data')

# Directory to save outputs.
output_dir = os.path.join(root_dir,'results','BodyPart_GazeTime_Tobii')


# Information and settings about densepose files
# Pickle file that contains denspose results. 
pickle_file = os.path.join(features_dir,'%s_densepose.pkl')

# denspose prediction threshold.
pred_score_thrs = 0.8

# Videos to compute gaze times.
vidclips = ['9.mp4', '12.mp4', '13.mp4', '14.mp4', '19.mp4', '21.mp4', '27.mp4','28.mp4','101.mp4', '102.mp4', '103.mp4', '104.mp4', '105.mp4', '106.mp4', '107.mp4','108.mp4', '109.mp4', '110.mp4']
gaze_radius = 15.5 # Half of 1 degree visual angle for Tobii

# set up the output directory
makedir(output_dir,sysexit=True)


for vii, vid_ii in enumerate(vidclips):
    
    print('Processing video file: %s'%vid_ii)

    # load downsampled gaze data, prepared in 02* scripts
    data_file = os.path.join(prepdata_dir,f'timebinned_data_{vid_ii}.pkl')
    vid_etdata = ETdata(data_file=data_file,stim_dir=stim_dir)
    ngroups = vid_etdata.data_ngroups


    # ----- Load some information about the video -----
    nframes = vid_etdata.stim_mediainfo['nframes']
    frame_width = vid_etdata.stim_mediainfo['frame_width']
    frame_height = vid_etdata.stim_mediainfo['frame_height']
        
    # ----- Load densepose results -----
    vid_basename = os.path.splitext(vid_etdata.stim_videoname)[0]
    with open(pickle_file%vid_basename,'rb') as pf:
        densepose_results = pickle.load(pf)
    
    et_xy = np.concatenate(vid_etdata.data, axis=0)
    subjs = np.hstack(vid_etdata.data_subjs)
    subj_groups = np.hstack(vid_etdata.data_subjs_group)
    
    gaze2screen = np.zeros(et_xy.shape[:2], dtype=bool)
    gaze2head = np.zeros(et_xy.shape[:2], dtype=bool)
    gaze2hands = np.zeros(et_xy.shape[:2], dtype=bool)
    gaze2wbody = np.zeros(et_xy.shape[:2], dtype=bool) # whole body 
    gaze2nonbody = np.zeros(et_xy.shape[:2], dtype=bool) # non-body = background
    gaze2nonheadbody = np.zeros(et_xy.shape[:2], dtype=bool)

    for fii in tqdm(range(nframes)):

        frame_name = f'frame_{fii+1}'
        frame_results = densepose_results.get(frame_name, None)
        
        if frame_results is not None:
            img_layers, _ = densepose_results2bodypart_layers(frame_results,[frame_height,frame_width],pred_score_thrs=pred_score_thrs)
            # if pred_score_thrs is larger than threshold used in densepose feature extraction then img_layers might return None. 
            if img_layers is not None: 
                # head:3, hands:2, other-body-parts:1, background:0
                img_layers_red = densepose_facehand_layers(img_layers) #.astype(float)
                img_layers_nonheadbody = np.logical_or(img_layers_red==1,img_layers_red==2)
                        
        
        for sii, subj_ii in enumerate(subjs):
            if fii >= et_xy[sii].shape[0]: # recording is shorter than video duration.
                continue

            this_et_xy = et_xy[sii,fii,:]
            
            if np.isnan(np.sum(this_et_xy)):
                continue
            if this_et_xy[0] < 0 or this_et_xy[1] < 0: 
                continue
            else:
                gaze2screen[sii,fii] = True # on-screen gaze. 
                gaze2nonbody[sii,fii] = True # assume for now. If not True, corrected below in gaze-to-any-part of the body condition. 
                if frame_results is not None and img_layers_red is not None:
                    
                    this_ET_mask = radmask(this_et_xy, gaze_radius, img_layers_red.shape)
                    
                    # gaze-to-head:
                    if np.any(img_layers_red==3):
                        overlap_TF = np.logical_and(img_layers_red==3, this_ET_mask).any()
                        if overlap_TF:
                            gaze2head[sii,fii] = True
                            
                    # gaze-to-hands:
                    if np.any(img_layers_red==2):
                        overlap_TF = np.logical_and(img_layers_red==2, this_ET_mask).any()
                        if overlap_TF:
                            gaze2hands[sii,fii] = True
                            
                    # gaze-to-any-part of the body:
                    if np.any(img_layers_red>0):
                        overlap_TF = np.logical_and(img_layers_red>0, this_ET_mask).any()
                        if overlap_TF:
                            gaze2wbody[sii,fii] = True
                            gaze2nonbody[sii,fii] = False 
                            
                    # gaze-to-any non-head part of the body:
                    if np.any(img_layers_red==3): # could reduce false positive body-part detections. 
                        overlap_TF = np.logical_and(img_layers_nonheadbody, this_ET_mask).any()
                        if overlap_TF:
                            gaze2nonheadbody[sii,fii] = True

                            
    # save this video outputs
    pkl_output_file = os.path.join(output_dir,f'{vid_ii}_fixtime_coeff_{gaze_radius}.pkl')
    data2out = {'onscreen':gaze2screen, 'heads':gaze2head, 'hands':gaze2hands, 'wbody':gaze2wbody, 
                'nonbody':gaze2nonbody, 'nonheadbody':gaze2nonheadbody, 
                'subj_info':subjs, 'subj_groups': subj_groups}
    with open(pkl_output_file,'wb') as fp:
        pickle.dump(data2out,fp)
    
    vidname_save = os.path.splitext(vid_ii)[0]
    savemat(os.path.join(output_dir,f'{vidname_save}_fixtime_coeff_{gaze_radius}.mat'), data2out)

    