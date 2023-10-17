#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

adapted from https://github.com/adolphslab/adolphslab_eyetracking
"""

import os
from alabeye.etdata import ETdata

# Main directory for experiment data.
root_dir = 'path-to-dir'

# HDF file of gaze data
hdf_file_list = [ 'ET_RemotePhone_session1.hdf5', 'ET_RemotePhone_session2.hdf5',
                  'ET_RemotePhone_session3.hdf5', 'ET_RemotePhone_session4.hdf5',
                  'ET_RemotePhone_session5.hdf5', 'ET_RemotePhone_session6.hdf5',
                  'ET_RemotePhone_session7.hdf5', 'ET_RemotePhone_session8.hdf5',
                  'ET_RemotePhone_session9.hdf5', 'ET_RemotePhone_session10.hdf5']


subj_info_file = os.path.join(root_dir, 'participants_info.csv')
stim_dir = 'path-to-stim'


for h5ii in hdf_file_list:
    
    sessionName = h5ii.split('_')[2]
    hdf_file = os.path.join(root_dir,'remotephone', h5ii)
    timebin_outdir = os.path.join(root_dir,'remotephone', 'down2frame_' + sessionName )
    
    etdata_init = ETdata(data_file=hdf_file)

    for task_ii in etdata_init.available_tasks:
        
        print(f'processing {task_ii}')
        etdata_task = ETdata(taskname=task_ii, data_file=hdf_file, 
                             subj_info_file=subj_info_file, use_subj_info=True,
                             stim_dir=stim_dir)


        # default is to bin to frame duration.
        etdata_task.get_timebinned_data(split_groups=True, fix_length=True,
                                         save_output=True, output_dir=timebin_outdir,
                                         rm_subj_withhalfmissingdata=True)
        
