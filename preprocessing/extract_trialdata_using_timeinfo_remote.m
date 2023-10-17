function data = extract_trialdata_using_timeinfo_remote(gazedir, session, block, trialtype)
%
% -----------------------------------------------------------------------
% This function takes gaze data from each block and trialtype 
% and re-organize by trials based on trial onset timestamps from metadata.json
%
% gazedir = directory path to psvr.mat
% block = block name as string (e.g. 'block1', 'block2')
% trialtype = trial type as string (e.g. 'dot', 'image_free_view')
%
% written by NK (Oct 2022)
% -----------------------------------------------------------------------

% read gaze data split by block
gazedata = load(fullfile(gazedir, 'psvr.mat'), block);
gazedata = gazedata.(block).(trialtype);

gazedata = table2cell(gazedata);

% extract time info from metadata
raw = fileread(fullfile(gazedir, ['autism_session' num2str(session) '_20211218_' block '_' trialtype '_0'], 'metadata.json'));
expression = '(?<=meta_data:)(.+?)(?=(meta_data:))';
[trials, lastTrial] = regexp(raw,expression,'tokens','split');
trials{end+1} = lastTrial(end);
t_onset = [];
for t=1:size(trials,2) %length(trials)
    raw_trial = trials{t}{1};
    expression = 'start_time_millis       : (\w+)';
    starttime = regexp(raw_trial,expression,'tokens');
    starttime = starttime{1}{1};
    t_onset = [t_onset; str2double(starttime)];
end

gazeidx = cell(size(trials,2),1);
for t=1:size(trials,2)
    if t==size(trials,2)
        gazeidx{t} = find([gazedata{:,3}] >= t_onset(t));
    else
        gazeidx{t} = find([gazedata{:,3}] >= t_onset(t) & [gazedata{:,3}] < t_onset(t+1));
    end
end


data = cell(size(trials,2),1);
for t=1:size(trials,2)     
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
