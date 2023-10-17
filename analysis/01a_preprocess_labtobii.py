#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Load data mat files and re-store in HDF file format.  
This version saves filtered data based on calibration quality in addition to non-filtered 

adapted and modified from https://github.com/adolphslab/adolphslab_eyetracking

"""


import os
import numpy as np
import pandas as pd
import pims

import json

from alabeye.io.loaddata import loadmat
from alabeye.stats import nanmasked_mean


# Directory of collected eye-tracking data (location of Tobii ET data)
exp_folder = 'path-to-data'

# Directory of stimulus videos
stim_dir = 'path-to-stim'

# Output file to save processed ET data. 
hdf_output_file = os.path.join(exp_folder, 'tobii', 'ET_Tobii.hdf5')
hdf_output_dir = os.path.dirname(hdf_output_file)

# saving two versions of output with and without data exclusion based on calibration quality
clean_hdf_output_file = os.path.join(exp_folder, 'tobii', 'ET_Tobii_clean.hdf5')
clean_hdf_output_dir = os.path.dirname(clean_hdf_output_file)

# if ET data contains very large number of subjects or very long tracking durations, 
# it might be better to compress hdf file. Compressed size is nearly the half of non-compressed size. 
# But compressing/uncompressing is a slow process.
compress_hdf = False

# Walk through subject-wise data folders. 
subj_files = []
for filename in os.listdir(exp_folder):
    if os.path.isdir(os.path.join(exp_folder,filename)):
        subj_files.append(filename)
subj_files = sorted(subj_files)

first_entry = True # to open an HDF file. 
first_entry_clean = True 
video_metainfo = {}

for fname_ii in subj_files:
    
    subj = fname_ii.split('_')[0]
    print(subj)
    if len(subj) != 6 or ' ' in subj:
        raise ValueError(f'Unexpected subj ID format: {subj}')
    
    # Read .mat file to load experiment data. 
    fname_ii = '%s_allmovies_gaze' % subj
    gaze_matfile = loadmat(os.path.join(exp_folder,subj,fname_ii))

    # Load raw gaze data for all movies. 
    gaze_data_all = gaze_matfile['AllData']
    vids_this = list(gaze_data_all.keys())
    
    # Load experiment parameters to rescale raw gaze data to movie frame size.
    params = gaze_matfile['Params']
    screenrect = params['ScreenRect']
    placeholder = params['placeholder']
    
    # from saved screen dimensions (experiment raw files)
    movHshown = placeholder[3] - placeholder[1]
    movWshown = placeholder[2] - placeholder[0]
    scale_factor = screenrect[3]/movHshown
    sub_factor = (screenrect[3]/movHshown - 1)/2
    
    
    # ---------------------------------------------------------------------
    # calibration metric
    # ---------------------------------------------------------------------
    screenwidth = screenrect[2]
    screenheight = screenrect[3]
    mon_width = 52.7 # tobii monitor width = 52.7cm
    view_distance_default = 62.5 # Tobii viewing distance = 60-65 cm
    calib_data = gaze_matfile['CalibrationResults']
    
    # loop through all videos
    allvids_caliberror_px_bypoint = pd.DataFrame()
    allvids_caliberror_ppd_bypoint = pd.DataFrame()
    allvids_lefteye_caliberror_px  = []
    allvids_righteye_caliberror_px = []
    allvids_botheyes_caliberror_px = []
    allvids_lefteye_caliberror_ppd  = []
    allvids_righteye_caliberror_ppd = []
    allvids_botheyes_caliberror_ppd = []    
    
    for vid_ii in vids_this:
        if type(params[vid_ii]['distance']) == float:
            params[vid_ii]['distance'] = np.array([params[vid_ii]['distance']])
            
        if len(params[vid_ii]['distance']) == 0:
            view_distance = view_distance_default           
        else:
            view_distance = np.mean(params[vid_ii]['distance'])    
        
        # pixel per degree
        ppd = screenwidth / (2* np.arctan(mon_width/(2*view_distance))*180/ np.pi);
        
        # extract calibration data and calculate deviations 
        if len(calib_data[vid_ii]['points_lefteye']) == 0:
            pass
        else:
            points = calib_data[vid_ii]['points_lefteye'].keys()
            assert len(calib_data[vid_ii]['points_lefteye']) == len(calib_data[vid_ii]['points_righteye'])
            points = calib_data[vid_ii]['points_lefteye'].keys()
        
            true_locations = calib_data[vid_ii]['points_on_display']
            calib_lefteye = calib_data[vid_ii]['points_lefteye']
            calib_righteye = calib_data[vid_ii]['points_righteye']
        
            lefteye_allpoints = []
            righteye_allpoints = []
            botheyes_allpoints = []
        
            for pii, pointname in enumerate(points):
                true_x = true_locations[pii][0] * screenwidth 
                true_y = true_locations[pii][1] * screenheight
                
                lefteye = calib_lefteye[pointname][11:,]
                righteye = calib_righteye[pointname][11:,]
                botheyes_x = nanmasked_mean([lefteye[:,0], righteye[:,0]], axis=0)
                botheyes_y = nanmasked_mean([lefteye[:,1], righteye[:,1]], axis=0)
        
                lefteye_dist  = np.sqrt((lefteye[:,0] * screenwidth - true_x)**2 + (lefteye[:,1] * screenheight - true_y)**2)
                righteye_dist = np.sqrt((righteye[:,0] * screenwidth - true_x)**2 + (righteye[:,1] * screenheight - true_y)**2)
                botheyes_dist = np.sqrt((botheyes_x * screenwidth - true_x)**2 + (botheyes_y * screenheight - true_y)**2)
        
                lefteye_allpoints.append(np.mean(lefteye_dist))
                righteye_allpoints.append(np.mean(righteye_dist))
                botheyes_allpoints.append(np.mean(botheyes_dist))
        
            allvids_lefteye_caliberror_px.append(np.mean(lefteye_allpoints))
            allvids_righteye_caliberror_px.append(np.mean(righteye_allpoints))
            allvids_botheyes_caliberror_px.append(np.mean(botheyes_allpoints))
            allvids_lefteye_caliberror_ppd.append(np.mean(lefteye_allpoints) / ppd)
            allvids_righteye_caliberror_ppd.append(np.mean(righteye_allpoints) / ppd)
            allvids_botheyes_caliberror_ppd.append(np.mean(botheyes_allpoints) / ppd)
        
            allvids_caliberror_px_bypoint[vid_ii] = [lefteye_allpoints, righteye_allpoints, botheyes_allpoints]
            allvids_caliberror_ppd_bypoint[vid_ii] = [lefteye_allpoints, righteye_allpoints, botheyes_allpoints]   
            # by-point data currently not being saved                                  


    to_compare = [np.mean(allvids_lefteye_caliberror_px), np.mean(allvids_righteye_caliberror_px), np.mean(allvids_botheyes_caliberror_px)]
    to_compare_TF = to_compare == np.min(to_compare) # 0 = left, 1 = right, 2 = both 
    
    # ---------------------------------------------------------------------
    # collect rescaled gaze data
    # ---------------------------------------------------------------------        

    for vid_ii in vids_this:

        # Rescale raw gaze data to movie frame size and save to an hdf file. 
        gaze_raw = gaze_data_all[vid_ii]
        
        data_len = len(gaze_raw['device_timestamp'])
        extra_info = {}
        for key,val in gaze_raw.items():
            if isinstance(val,int):
                extra_info[key] = val
            elif len(val) != data_len:
                extra_info[key] = val
        
        # delete extra_info from gaze_raw before converting to dataframe.
        for key in extra_info.keys():
            gaze_raw.pop(key)
        
        # Actual file name of the video.
        vid_fname = f"{vid_ii.replace('mv','')}.mp4"
        if vid_fname not in video_metainfo:
            # define VideoReader object and get metadata from the video.
            video_file = os.path.join(stim_dir,vid_fname)
            if not os.path.isfile(video_file):
                raise SystemExit(f'Could not find the video file:\n{video_file}')
            
            vr = pims.PyAVReaderTimed(video_file)
            frame_width, frame_height = vr.frame_shape[1], vr.frame_shape[0]
            vid_duration, nframes, vid_fps = vr.duration, len(vr), vr.frame_rate

            video_metainfo[vid_fname] = {'frame_width':frame_width, 
                                         'frame_height':frame_height,
                                         'vid_duration':vid_duration,
                                         'nframes':nframes,'vid_fps':vid_fps}

        # take and save relevant info.
        gaze_data_df = pd.DataFrame.from_dict(gaze_raw,dtype=float)
        
        frame_width = video_metainfo[vid_fname]['frame_width']
        frame_height = video_metainfo[vid_fname]['frame_height']
        
        # the same as first averaging left and right eyes and scaling to frame size. 
        gaze_data_df[ ['left_x_coord', 'right_x_coord'] ] = ( (gaze_data_df[ ['left_x_coord', 'right_x_coord'] ]*scale_factor - sub_factor)*frame_width)
        gaze_data_df[ ['left_y_coord', 'right_y_coord'] ] = ( (gaze_data_df[ ['left_y_coord', 'right_y_coord'] ]*scale_factor - sub_factor)*frame_height)
        
        if to_compare_TF[0] == True:
            gaze_x = gaze_data_df['left_x_coord'].values
            gaze_y = gaze_data_df['left_y_coord'].values
            summary_caliberror_ppd = to_compare[0] / ppd
            summary_caliberror_percent = to_compare[0] / movHshown
        elif to_compare_TF[1] == True:
            gaze_x = gaze_data_df['right_x_coord'].values
            gaze_y = gaze_data_df['right_y_coord'].values
            summary_caliberror_ppd = to_compare[1] / ppd
            summary_caliberror_percent = to_compare[1] / movHshown
        elif to_compare_TF[2] == True:
            gaze_x = nanmasked_mean(gaze_data_df[ ['left_x_coord', 'right_x_coord'] ].values, axis=1)
            gaze_y = nanmasked_mean(gaze_data_df[ ['left_y_coord', 'right_y_coord'] ].values, axis=1)
            summary_caliberror_ppd = to_compare[2] / ppd
            summary_caliberror_percent = to_compare[2] / movHshown
        
        # it is tricky to average pupil diameters from two eyes. Let's keep both for now. 
        # pupil_diameter = nanmasked_mean(gaze_data_df[ ['left_pupil_diameter', 'right_pupil_diameter'] ].values, axis=1)

        # put data into a format that we used in previous studies. 
        ETdata_df = pd.DataFrame({'RecTime': gaze_data_df['device_timestamp'].values, 
                                  'GazeX': gaze_x, 'GazeY': gaze_y, 
                                  'PupilDiameter_left':gaze_data_df['left_pupil_diameter'].values,
                                  'PupilDiameter_right':gaze_data_df['right_pupil_diameter'].values,
                                  'GazeErrorPpd': summary_caliberror_ppd, 
                                  'GazeErrorPercent': summary_caliberror_percent})
        
        # changing mv to .mp4 to avoid confusion 
        vid_ii = vid_fname
        
        session_ID = '%s/%s'%(subj,vid_ii) # using original stimulus name used in .mat files, not vid_fname.
        if first_entry:
            if compress_hdf:
                ETdata_df.to_hdf(hdf_output_file,session_ID,mode='w',complevel=9, complib='zlib')
            else:    
                ETdata_df.to_hdf(hdf_output_file,session_ID,mode='w')
            first_entry=False
        else:
            if compress_hdf:
                ETdata_df.to_hdf(hdf_output_file,session_ID,mode='a',complevel=9, complib='zlib')
            else:
                ETdata_df.to_hdf(hdf_output_file,session_ID,mode='a')
        
        # include in the 'clean' dataset only if calibration error < 2 deg
        if summary_caliberror_ppd < 2:
            if first_entry_clean:
                if compress_hdf:
                    ETdata_df.to_hdf(clean_hdf_output_file,session_ID,mode='w',complevel=9, complib='zlib')
                else:    
                    ETdata_df.to_hdf(clean_hdf_output_file,session_ID,mode='w')
                first_entry_clean=False
            else:
                if compress_hdf:
                    ETdata_df.to_hdf(clean_hdf_output_file,session_ID,mode='a',complevel=9, complib='zlib')
                else:
                    ETdata_df.to_hdf(clean_hdf_output_file,session_ID,mode='a')

        # save video meta info
        with open('stimname_mvid_ap.json', 'w') as f:
             json.dump(video_metainfo, f)
        

