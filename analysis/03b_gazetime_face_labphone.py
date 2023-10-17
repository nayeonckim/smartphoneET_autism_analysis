#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute gaze time to face area and center.
adapted and modified from https://github.com/adolphslab/adolphslab_eyetracking

"""
import sys
import os
import numpy as np
import pickle
from tqdm import tqdm

# The Adolphs Lab data analysis scripts for eye tracking are provided under alabeye package  
from alabeye.etdata import ETdata, makedir
from alabeye.io.loadstimfeatures import get_faceareas_simple, radmask, sort_facedetect_by_salience, get_faceareas

from scipy.spatial.distance import cdist

# if matlab version of output is needed
from scipy.io import savemat 

def point_in_rectangle(x, y, rect_x, rect_y, rect_width, rect_height):
    if x >= rect_x and x <= rect_x + rect_width:
        if y >= rect_y and y <= rect_y + rect_height:
            return True
    return False


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
prepdata_dir = os.path.join(root_dir,'labphone', 'down2frame_data')

# Directory to save outputs
output_dir = os.path.join(root_dir,'results','Face_GazeTime_LabPhone')

# Information and settings about face detection files
# Pickle file that contains retinaface results. 
pickle_file_facedetect = os.path.join(features_dir,'%s_facedetect.pkl')

# Videos to compute gaze times.
vidclips = ['9.mp4', '12.mp4', '13.mp4', '14.mp4', '19.mp4', '21.mp4', '27.mp4','28.mp4','101.mp4', '102.mp4', '103.mp4', '104.mp4', '105.mp4', '106.mp4', '107.mp4','108.mp4', '109.mp4', '110.mp4']
gaze_radius = 28.5 # Half of 1 degree visual angle for phone

# define center area to calculate center bias 
center_x_min = 506
center_y_min = 226
center_x_len = 268
center_y_len = 268

# set up the output directory
makedir(output_dir,sysexit=False)


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


    # ----- Load facedetect / retinaface results -----
    vid_basename = os.path.splitext(vid_etdata.stim_videoname)[0]
    with open(pickle_file_facedetect%vid_basename,'rb') as pf:
        facedetect_results = pickle.load(pf)

    et_xy = np.concatenate(vid_etdata.data, axis=0)
    subj_lists = vid_etdata.data_subjs
    subj_groups = np.hstack(vid_etdata.data_subjs_group)
    subjs = np.hstack(subj_lists)

    
    face_use_info = []
    face_size_info = []
    
    face_gaze = np.zeros(et_xy.shape[:2],dtype=np.int8)
    face_gaze_dists = np.full((*et_xy.shape[:2],3),np.nan) # 3 distances. 
    face_gaze_dists_valid = np.zeros(et_xy.shape[:2],dtype=bool)
    face_num = np.zeros(et_xy.shape[:2],dtype=int)
    
    center_gaze = np.zeros(et_xy.shape[:2],dtype=np.int8)

    for fii in tqdm(range(nframes)):

        # --- Load ET data for both groups ---    
        et_bin = et_xy[:,fii,:]
              
        # --- center gaze (gaze disc is not used) --- 
        for sub_ii, sub_et_ii in enumerate(et_bin):
            if sub_et_ii[0] > 0 and sub_et_ii[1] > 0 and not np.isnan(sub_et_ii.mean()): 
               if point_in_rectangle(sub_et_ii[0], sub_et_ii[1], center_x_min, center_y_min, center_x_len, center_y_len):
                  center_gaze[sub_ii, fii] = 1  
            else:
               center_gaze[sub_ii, fii] = -1  
        
        
        # --- center gaze (gaze disc is used) ---  
        scan_area = []
        for sub_ii, sub_et_ii in enumerate(et_bin):
            if sub_et_ii[0] > 0 and sub_et_ii[1] > 0 and not np.isnan(sub_et_ii.mean()):               
                this_scan_mask = radmask(sub_et_ii, gaze_radius, [frame_height,frame_width])
                scan_area.append(this_scan_mask)
            else:
                scan_area.append(None)

        # --- Load face-detection results ---
        frame_name = f'frame_{fii+1}'
        frame_results_facedetect = facedetect_results.get(frame_name,None)
        
        face_areas = None
        face_use = False
        face_size = False
        if frame_results_facedetect is not None:

            frame_results_facedetect = sort_facedetect_by_salience(frame_results_facedetect,scan_area,frame_height,frame_width)
            face_areas, landmarks = get_faceareas_simple(frame_results_facedetect,frame_height,frame_width,detection_thrs=0.5)

            # assess face size to determine whether a face size is large enough to reliably measure distances from landmarks. 
            diff_faces = np.unique(face_areas[face_areas>0])
            if len(diff_faces) != len(landmarks):
                print('\n Overlapping face areas!')
                print(diff_faces)
                landmarks = [ landmarks[int(of_ii)] for of_ii in diff_faces-1 ]
                diff_faces = np.arange(1,len(landmarks)+1).astype(float)
                print(diff_faces)
            
            assert len(diff_faces) == len(landmarks)

            face_use = np.ones((len(diff_faces))).astype(bool)
            for la_ii, landmark in enumerate(landmarks):
                land_dists = cdist( landmark[[0,1],:], landmark[[2,3],:] ) # distance between eyes and nose or mouth. 

                if (land_dists < gaze_radius*2.).any():
                    face_use[la_ii] = False

            face_size = np.zeros((len(diff_faces)))
            for face_cnt,face_ii in enumerate(diff_faces):
                this_face = face_areas==face_ii
                face_size[face_cnt] = np.sum(this_face>0)

        # keep face use info.
        face_use_info.append(face_use)
        face_size_info.append(face_size)

        
        # Count face gaze and measure gaze to landmark distances.
        if not face_areas is None:
            for sc_cnt,sc_ii in enumerate(scan_area):
                if sc_ii is None:
                    face_gaze[sc_cnt,fii] = -1

                else:
                    for face_cnt,face_ii in enumerate(diff_faces):
                        this_face = face_areas==face_ii

                        this_landmarks = landmarks[face_cnt]
                        overlap_TF = np.logical_and(this_face, sc_ii).any()

                        if overlap_TF:
                            face_gaze[sc_cnt,fii] = 1
        
                            dist_dum = cdist([et_bin[sc_cnt,:]],this_landmarks).squeeze()
                            face_gaze_dists[sc_cnt,fii,:] = [ dist_dum[[0,1]].min(), dist_dum[2], dist_dum[3] ]
                            face_gaze_dists_valid[sc_cnt,fii] = face_use[int(face_ii)-1]
                            assert face_cnt == int(face_ii)-1, 'Expected pattern...'
                            face_num[sc_cnt,fii] = face_ii
                            break


        else: # onscreen fixation when there is no face on the screen. 
            for sc_cnt,sc_ii in enumerate(scan_area):
                if sc_ii is None:
                    face_gaze[sc_cnt,fii] = -1


    pkl_output_file = os.path.join(output_dir,f'{vid_ii}_fixtime_coeff_{gaze_radius}.pkl')
    data2out = {'face_use_info':face_use_info, 'face_size_info':face_size_info, 
                'subjs':subjs, 'subj_lists':subj_lists, 'subj_groups':subj_groups,
                'face_gaze':face_gaze, 'center_gaze':center_gaze, 
                'face_gaze_dists':face_gaze_dists, 'face_gaze_dists_valid':face_gaze_dists_valid,
                'face_num':face_num}
    with open(pkl_output_file,'wb') as fp:
        pickle.dump(data2out,fp)

    savemat(os.path.join(output_dir,f'{vid_basename}_fixtime_coeff_{gaze_radius}.mat'), data2out)    
    