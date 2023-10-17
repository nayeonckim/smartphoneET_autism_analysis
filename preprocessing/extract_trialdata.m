function data = extract_trialdata(gazedir, block, trialtype)
%
% -----------------------------------------------------------------------
% This function takes gaze data from each block and trialtype 
% and re-organize by trials based on timing info recorded in psvr file
% *** Recommend using "extract_trialdata_using_timinginfo" isntead of this ***

% gazedir = directory path to psvr.mat
% block = block name as string (e.g. 'block1', 'block2')
% trialtype = trial type as string (e.g. 'dot', 'image_free_view')
%
% written by NK  (May 2022)
% -----------------------------------------------------------------------


% read gaze data split by block
gazedata = load(fullfile(gazedir, 'psvr.mat'), block);
gazedata = gazedata.(block).(trialtype);

gazedata = table2cell(gazedata);

% split individual trials
split_idx = find(diff([gazedata{:,4}]) < 0);
gazeidx = {};
numTrials = length(split_idx)+1;
for t=1:numTrials
    if t==1
        gazeidx{t} = 1:split_idx(t);
    elseif t==numTrials
        gazeidx{t} = (split_idx(end)+1):length(gazedata);
    else
        gazeidx{t} = (split_idx(t-1)+1):split_idx(t);
    end
end

data = cell(numTrials,1);
for t=1:numTrials     
    xx = cell2mat(gazedata(gazeidx{t},5));
    yy = cell2mat(gazedata(gazeidx{t},6));
    tt_orig = cell2mat(gazedata(gazeidx{t},3));  
    tt_since_onset = cell2mat(gazedata(gazeidx{t},4)); 
    
    data{t} = [xx, yy, tt_orig, tt_since_onset];
    
    if strcmp(trialtype,'dot')
        true_xx = cell2mat(gazedata(gazeidx{t},7));
        true_yy = cell2mat(gazedata(gazeidx{t},8));
        data{t} = [data{t}, true_xx, true_yy];
    end
    
end 
