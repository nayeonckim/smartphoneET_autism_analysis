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
hdf_file_list = [ 'ET_Tobii.hdf5' ]

subj_info_file = os.path.join(root_dir,'participants_info.csv')
stim_dir = 'path-to-stim'

timebin_outdir = os.path.join(root_dir,'tobii','down2frame_data')


for h5ii in hdf_file_list:
    
    hdf_file = os.path.join(root_dir,'tobii',h5ii)
    
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
        
