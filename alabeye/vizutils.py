#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""



import numpy as np
from PIL import Image


def rgb2gray3ch(rgb):
    gray1ch = np.dot(rgb[...,:3], [0.299, 0.587, 0.114])
    w, h = gray1ch.shape
    ret = np.empty((w, h, 3), dtype=np.uint8)
    ret[:, :, 2] =  ret[:, :, 1] =  ret[:, :, 0] =  gray1ch
    return ret


def convertMillis(millis):
    seconds=np.int((millis/1000)%60)
    minutes=np.int((millis/(1000*60))%60)
    
    return "{:02d}m {:02d}s".format(minutes, seconds)
    


def resize_image(im, size, filter=Image.ANTIALIAS):
    '''Resize an image and return its array representation.
    size=(width, height) or downsample rate. 

    Parameters
    ----------
    im : str, np.ndarray(uint8), or PIL.Image object
        The path to the image, an image array, or a loaded PIL.Image.
    size : tuple, (width, height) or scalar, downsample rate.
        The desired output image size.

    Returns
    -------
    arr : uint8 narray, (height, width, 3)
        The resized image array
    '''
    if isinstance(im, str):
        im = Image.open(im)
    elif isinstance(im, np.ndarray):
        im = Image.fromarray(im)

    if  isinstance(size, (int, float)):
        size_wh = (im.width // size, im.height // size)
    else:
        size_wh = size

    im_resized = im.resize(size_wh, filter)
    im_resized = np.asarray(im_resized)
    
    return im_resized



def radMask(center_xy,radius,array_shape):
    
    if np.array(center_xy).shape != (2,): 
        raise SystemExit('Problem in the input of radMask')
    
    col_val = center_xy[0] 
    row_val = center_xy[1] 
    n_rows,n_cols = array_shape
    rows,cols = np.ogrid[-row_val:n_rows-row_val,-col_val:n_cols-col_val]
    mask = cols*cols + rows*rows <= radius*radius
    return mask    


def get_boundingbox(arr_):
    rows = np.any(arr_, axis=1)
    cols = np.any(arr_, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    # [xmin, ymin, xmax, ymax] corresponding to [ left, top, right, bottom ]
    return np.array([cmin, rmin, cmax, rmax]) 



# used in fixation analyses for RetinaFace detections. 
def sort_facedetect_by_salience(this_frame_results,scan_area_Ctrl,frame_height,frame_width):
    n_subjs = []
    for hii in this_frame_results:
        face_areas = np.zeros((frame_height,frame_width)).astype(bool)

        b = np.maximum(hii.astype(int), 0)
        face_areas[b[1]:b[3],b[0]:b[2]] = True
        
        n_subjs.append( sum( np.logical_and(sii,face_areas).any() for sii in scan_area_Ctrl if (not sii is None) ) )
        
    this_frame_results = this_frame_results[np.argsort(n_subjs)[::-1]] # most salient to less. 
    return this_frame_results


def get_faceareas_simple(this_frame_results, frame_height, frame_width,
                         detection_thrs=0.5, sort_by_area=False,
                         return_corners=False, return_landmarks_full=False):

    # sort faces based on face size, rather than prediction probability.
    if sort_by_area:
        areas_ = []
        for hii in this_frame_results:
            areas_.append((hii[3]-hii[1])*(hii[2]-hii[0]))
        this_frame_results = this_frame_results[np.argsort(areas_)[::-1]]
    
    face_areas = np.zeros((frame_height,frame_width))

    if return_corners:
        corners = []
    
    face_num = 1
    landmarks_all = []
    for hii in this_frame_results:
        if hii[4] > detection_thrs:
            
            b = np.maximum(hii.astype(int), 0)
            if b[2]>frame_width: b[2] = frame_width-1 
            if b[3]>frame_height: b[3] = frame_height-1
            
            face_areas[b[1]:b[3],b[0]:b[2]] = face_num
            face_num += 1
            
            landmarks = b[5:].reshape(-1,2)
            #from viewers perspective: b5,b6: left eye; b7,b8: right eye;  
            # b9,b10: nose; b11,b12: left-side mouth; b13,b14: right-side mouth
            
            # control for out of frame detections. 
            landmarks[:, 0] = np.minimum(landmarks[:, 0], frame_width - 1)
            landmarks[:, 1] = np.minimum(landmarks[:, 1], frame_height - 1)
            
            assert landmarks.shape[0]==5
            # center_of_eyes = np.mean(landmarks[[0,1],:],axis=0).round().astype(int)
            mouth_center = np.mean(landmarks[[3,4],:],axis=0).round().astype(int)
            
            if return_landmarks_full:
                landmarks_all.append(np.vstack((landmarks)))
            else:
                landmarks_all.append(np.vstack((landmarks[:-2,:], mouth_center)))

            if return_corners:
                corners.append(b[0:4])
    
    if return_corners:
        return face_areas, landmarks_all, corners
    
    return face_areas, landmarks_all



# used in fixation analyses for RetinaFace detections. 
def get_faceareas(this_frame_results, frame_height, frame_width, detection_thrs=0.5):

    face_areas = np.zeros((frame_height,frame_width),dtype=bool)
    face_areas_eyes = np.zeros((frame_height,frame_width),dtype=bool)
    face_areas_mouth = np.zeros((frame_height,frame_width),dtype=bool)
    
    for hii in this_frame_results:
        if hii[4] > detection_thrs:
            
            b = np.maximum(hii.astype(int), 0)
            
            face_areas[b[1]:b[3],b[0]:b[2]] = True
            
            landmarks = b[5:].reshape(-1,2)
            #from viewers perspective: b5,b6: left eye; b7,b8: right eye;  
            # b9,b10: nose; b11,b12: left-side mouth; b13,b14: right-side mouth
            
            # control for out of frame detections. 
            landmarks[:, 0] = np.minimum(landmarks[:, 0], frame_width - 1)
            landmarks[:, 1] = np.minimum(landmarks[:, 1], frame_height - 1)
            
            box2 = [ np.min(landmarks[:,0]), np.min(landmarks[:,1]), 
                     np.max(landmarks[:,0]), np.max(landmarks[:,1]) ] 

            med_box = np.vstack((b[:4],box2)).mean(0).astype(int).reshape(-1,2)
            
            # Control for out of frame detections.
            med_box[:, 0] = np.minimum(med_box[:, 0], frame_width - 1)
            med_box[:, 1] = np.minimum(med_box[:, 1], frame_height - 1)
            med_box = med_box.flatten()

            # --- Following 2 parts will fail for non-vertical faces, 
            # but a good approximation for the Office videos. ---
            # upper - eye - 
            box_11 = [np.min([ med_box[0], med_box[2], landmarks[2,0] ]), np.min([med_box[1], landmarks[2,1]]), 
                      np.max([ med_box[0], med_box[2], landmarks[2,0] ]), np.max([med_box[1], landmarks[2,1]]) ]
            box_11 = np.maximum(box_11, 0)
            face_areas_eyes[box_11[1]:box_11[3],box_11[0]:box_11[2]] = True

            # lower - mouth - 
            box_22 = [np.min([ med_box[0], med_box[2], landmarks[2,0] ]), np.min([med_box[3], landmarks[2,1]]), 
                      np.max([ med_box[0], med_box[2], landmarks[2,0] ]), np.max([med_box[3], landmarks[2,1]]) ]
            box_22 = np.maximum(box_22, 0)
            face_areas_mouth[box_22[1]:box_22[3],box_22[0]:box_22[2]] = True

    return face_areas, face_areas_eyes, face_areas_mouth



# used in fixation analyses for RetinaFace detections. 
def get_faceareas_boxes(this_frame_results, frame_height, frame_width, detection_thrs=0.5):

    face_boxes = []
    face_boxes_eyes = []
    face_boxes_mouth = []
    
    for hii in this_frame_results:
        if hii[4] > detection_thrs:
            
            b = np.maximum(hii.astype(int), 0)
            
            face_boxes.append( [b[1], b[3], b[0], b[2]] )
            
            landmarks = b[5:].reshape(-1,2)
            #from viewers perspective: b5,b6: left eye; b7,b8: right eye;  
            # b9,b10: nose; b11,b12: left-side mouth; b13,b14: right-side mouth
            
            # control for out of frame detections. 
            landmarks[:, 0] = np.minimum(landmarks[:, 0], frame_width - 1)
            landmarks[:, 1] = np.minimum(landmarks[:, 1], frame_height - 1)
            
            box2 = [ np.min(landmarks[:,0]), np.min(landmarks[:,1]), 
                     np.max(landmarks[:,0]), np.max(landmarks[:,1]) ] 

            med_box = np.vstack((b[:4],box2)).mean(0).astype(int).reshape(-1,2)
            
            # Control for out of frame detections.
            med_box[:, 0] = np.minimum(med_box[:, 0], frame_width - 1)
            med_box[:, 1] = np.minimum(med_box[:, 1], frame_height - 1)
            med_box = med_box.flatten()

            # --- Following 2 parts will fail for non-vertical faces, 
            # but a good approximation for the Office videos. ---
            # upper - eye - 
            box_11 = [np.min([ med_box[0], med_box[2], landmarks[2,0] ]), np.min([med_box[1], landmarks[2,1]]), 
                      np.max([ med_box[0], med_box[2], landmarks[2,0] ]), np.max([med_box[1], landmarks[2,1]]) ]
            box_11 = np.maximum(box_11, 0)
            face_boxes_eyes.append([box_11[1], box_11[3], box_11[0], box_11[2]] )
            
            # lower - mouth - 
            box_22 = [np.min([ med_box[0], med_box[2], landmarks[2,0] ]), np.min([med_box[3], landmarks[2,1]]), 
                      np.max([ med_box[0], med_box[2], landmarks[2,0] ]), np.max([med_box[3], landmarks[2,1]]) ]
            box_22 = np.maximum(box_22, 0)
            face_boxes_mouth.append([box_22[1], box_22[3], box_22[0], box_22[2]] )
            
    return face_boxes, face_boxes_eyes, face_boxes_mouth




