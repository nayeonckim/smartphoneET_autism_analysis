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

# Directory of collected eye-tracking data 
exp_folder = 'path-to-data'

# Directory of stimulus videos
stim_dir = 'path-to-stim'

# Output file to save processed ET data. 
hdf_output_file = os.path.join(exp_folder,'labphone', 'ET_LabPhone.hdf5')
hdf_output_dir = os.path.dirname(hdf_output_file)

# saving two versions of output with and without data exclusion based on calibration quality
clean_hdf_output_file = os.path.join(exp_folder,'labphone', 'ET_LabPhone_clean.hdf5')
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
    fname_ii = '%s_allmovies_gaze_labphone' % subj
    gaze_matfile = loadmat(os.path.join(exp_folder,subj,fname_ii))

    # Load raw gaze data for all movies. 
    gaze_data_all = gaze_matfile['AllData']
    vids_this = list(gaze_data_all.keys())
    
    # Load experiment parameters to rescale raw gaze data to movie frame size.
    params = gaze_matfile['Params']
    screenrect = params['ScreenRect']
    placeholder = params['placeholder']
    
    # from saved screen dimensions 
    movHshown = placeholder[3] - placeholder[1]
    movWshown = placeholder[2] - placeholder[0]
    sub_factor = 120
    
    # ---------------------------------------------------------------------
    # calibration metric
    # ---------------------------------------------------------------------   

    screenwidth = screenrect[2]
    screenheight = screenrect[3]
    mon_width = 13.6 # phone monitor width = 13.6 cm 
    view_distance_default = 30.5 # phone viewing distance = 11-13 inches = 30.5 cm
    calib_data = gaze_matfile['CalibrationResults']
    
    # loop through all videos
    allvids_caliberror_px_bypoint = pd.DataFrame()
    allvids_caliberror_ppd_bypoint = pd.DataFrame()
    allvids_caliberror_px  = []
    allvids_caliberror_ppd  = []
    
    for vid_ii in vids_this:
    
        view_distance = view_distance_default
        # pixel per degree
        ppd = screenwidth / (2* np.arctan(mon_width/(2*view_distance))*180/ np.pi);
    
        calib_allpoints = []
    
        for pii in range(len(calib_data[vid_ii])):
            this_calib_data = pd.DataFrame(calib_data[vid_ii][pii], columns = ['calib_x', 'calib_y', 'device_timestamp', 'timestamp', 'true_x', 'true_y'])
            true_x = this_calib_data['true_x'][0]
            true_y = this_calib_data['true_y'][0]
            calib_x = np.array(this_calib_data[this_calib_data['timestamp'] < 1300]['calib_x']) 
            calib_y = np.array(this_calib_data[this_calib_data['timestamp'] < 1300]['calib_y']) 
    
            calib_dist = np.sqrt((calib_x - true_x)**2 + (calib_y - true_y)**2)
            calib_allpoints.append(np.mean(calib_dist))
    
        allvids_caliberror_px.append(np.mean(calib_allpoints))
        allvids_caliberror_ppd.append(np.mean(calib_allpoints) / ppd)
        
        if len(calib_allpoints) == 13:
            allvids_caliberror_px_bypoint[vid_ii] = calib_allpoints
            allvids_caliberror_ppd_bypoint[vid_ii] = calib_allpoints 
        else:
            print(subj + ' ' + vid_ii + ' has missing data')
    
    summary_caliberror_px = np.mean(allvids_caliberror_px)    
       
    # ---------------------------------------------------------------------
    # collect rescaled gaze data
    # ---------------------------------------------------------------------     
    
    for vid_ii in vids_this:

        # Rescale raw gaze data to movie frame size and save to an hdf file. 
        gaze_raw = gaze_data_all[vid_ii]
        
        data_len = len(gaze_raw['timestamp_res'])
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
        gaze_data_df['x_coord_res'] = ( (gaze_data_df['x_coord_res'] - sub_factor)/movWshown * frame_width)
        gaze_data_df['y_coord_res'] = ( gaze_data_df['y_coord_res'] / movHshown * frame_height)
        gaze_data_df['timestamp_res'] = ( gaze_data_df['timestamp_res'] - 0.2 )
        gaze_data_df = gaze_data_df[gaze_data_df['timestamp_res'] > 0]

        gaze_x = gaze_data_df['x_coord_res'].values
        gaze_y = gaze_data_df['y_coord_res'].values
        
        summary_caliberror_ppd = summary_caliberror_px / ppd
        summary_caliberror_percent = summary_caliberror_px / movHshown
        

        # put data into a format that we used in previous studies. 
        ETdata_df = pd.DataFrame({'RecTime': gaze_data_df['timestamp_res'].values, 
                                  'GazeX': gaze_x, 'GazeY': gaze_y,
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
        
        # include in the 'clean' dataset only if calibration error < 1.48 deg
        if summary_caliberror_ppd < 1.48:
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
        with open('stimname_map.json', 'w') as f:
             json.dump(video_metainfo, f)

        
