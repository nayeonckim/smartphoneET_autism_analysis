%
% -----------------------------------------------------------------------
% This script organizes eyedata and trial info from smartphone data 
% REMOTE PHONE + MOVIE-VIEW BLOCK (BLOCKS 2-4)
% 
% Requires: 
% - extract_trialdata.m (previous version) 
% - extract_trialdata_using_timeinfo_remote.m (more updated version)
% - metadata.json for each block and for each subj
% - psvr.mat for each subj 
%
% OUTPUT: a MAT file for each block containing 
% - calibdata: calibration data ("dot" trial)
% - eyedata: gaze data per trial
% - images: stimulus image filename per trial
%
% written by NK (Oct 2022)
% -----------------------------------------------------------------------

% Match Phone ID (automatically generated while saving to the server) to subject ID 
idmatch = readtable('idmatch_all.csv'); 
localdir = 'path-to-raw-data';

for sessionNum=1:10
    
    rawdata.dir       = fullfile(localdir, ['autism_session' num2str(sessionNum) '_20211218']);
    rawdata.IDs     = dir([rawdata.dir '/P_*']);
    
    % Match Phone ID (automatically generated while saving to the server) to subject ID 
    IDs = {};
    for s=1:length(rawdata.IDs)
        IDs{s,1} = rawdata.IDs(s).name;

        for tt=1:size(idmatch,1)
            if strcmp(rawdata.IDs(s).name, idmatch.PhoneID{tt})
                subjID = idmatch.SubjectID{tt};
            end
        end

        IDs{s,2} = subjID;
        clear subjID
    end

    %%

    summary_table_allsubj = {};
    for s=1:size(IDs,1)
        phoneID = IDs{s,1};
        subjID = IDs{s,2};

        % define directories
        data_dir_gaze = fullfile(rawdata.dir, phoneID);

        for c=1:3
            b = c+1; % block number = 2 to 4
            blockname = [ 'block' num2str(b) ];
            data_dir_meta = fullfile(data_dir_gaze, ['autism_session' num2str(sessionNum) '_20211218_' blockname '_video_free_view_0']);

            % organize metadata + gaze 
            if exist(data_dir_meta, 'dir')    

                % extract gaze data organized by trial 
                calibdata = extract_trialdata(data_dir_gaze, blockname, 'dot');
                eyedata   = extract_trialdata_using_timeinfo_remote(data_dir_gaze, sessionNum, blockname, 'video_free_view');

                % read metadata
                fname = fullfile(data_dir_meta, 'metadata.json');
                raw = fileread(fname);
                
                % extract trial info from metadata
                expression = '(?<=meta_data:)(.+?)(?=(meta_data:))';
                [trials, lastTrial] = regexp(raw,expression,'tokens','split');
                trials{end+1} = lastTrial(end);

                summary_table = cell(length(trials), 6); % phoneID, subjID, session, trial index, videofile, start_time
                images = cell(length(trials), 2);
                for t=1:length(trials)
                    raw_trial = trials{t}{1};
                    
                    % double check if participant ID matches
                    expression = 'participant_id: "(\w+)"';
                    phoneID_from_file = regexp(raw_trial,expression,'tokens');
                    if ~strcmp(phoneID_from_file{1}{1}, phoneID)
                        disp('WARNING: Phone ID does not match')
                    end

                    summary_table{t,1} = phoneID;
                    summary_table{t,2} = subjID;
                    summary_table{t,3} = sessionNum;

                    summary_table{t,4} = t;

                    expression = ['path  : "/storage/emulated/0/Android/media/autism_longitudinal_data/session' num2str(sessionNum) '/videos/block' num2str(c) '/(\w+).mp4'];
                    imagename = regexp(raw_trial,expression,'tokens');
                    summary_table{t,5} = [imagename{1}{1} '.mp4'];

                    images{t,1} = [imagename{1}{1} '.mp4'];

                    expression = 'start_time_millis       : (\w+)';
                    starttime = regexp(raw_trial,expression,'tokens');
                    starttime = starttime{1}{1};
                    summary_table{t,6} = starttime;

                    if isempty(eyedata{t})
                        tt = [];
                        disp(['WARNING for ' subjID ' (' phoneID ') ' blockname ' trial '  num2str(t) ' --CHECK'])
                    else
                        tt = eyedata{t}(:,3) - str2double(starttime);
                    end

                    eyedata{t} = [eyedata{t}, tt]; 


                    clear imagename starttime                 
                end

                if ~exist(fullfile(data_dir_gaze, [ blockname '_video_free_view.mat']), 'file')
                    fprintf(['saving .mat file for ' phoneID ' '  subjID '... \n']);
                    save(fullfile(data_dir_gaze, [ blockname '_video_free_view.mat']), 'calibdata', 'eyedata', 'images');
                    T = cell2table(summary_table, "VariableNames", ["phone_id" "subj_id" "session_num" "trial_num" "stim_file" "start_time"]);
                    writetable(T, fullfile(data_dir_gaze,  [ blockname '_summary.csv']))  
                end

                summary_table_allsubj = [summary_table_allsubj; summary_table];
            end  
        end
    end

    fprintf(['saving summary file for session ' num2str(sessionNum) '... \n']);
    T = cell2table(summary_table_allsubj, "VariableNames", ["phone_id" "subj_id" "session_num" "trial_num" "stim_file" "start_time"]);
    writetable(T, fullfile(localdir,  ['allsubj_video_session' num2str(sessionNum) '_summary.csv'])) 

end