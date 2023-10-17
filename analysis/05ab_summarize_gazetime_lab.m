% this script reorganizes the results files to generate figures 


rootdir = 'results';

videos = {'9'; '12'; '13'; '14'; '19'; '21'; '27'; '28';'101'; '102'; '103'; '104'; '105'; '106'; '107'; '108'; '109'; '110'};
bodyfeatures = {'wbody'; 'nonbody'};
facefeatures = {'face_gaze'; 'center_gaze'};


%% index video frames without any human body parts 

% assumption: when no one has gaze to wbody features across the two methods, it means there's no body feature on the frame. 

bodyparts_collect = [];

for vd=1:numel(videos)
    vdnum = videos{vd};
    vdname = ['mv' vdnum];
    data1 = load(fullfile(rootdir, 'BodyPart_GazeTime_Tobii', [vdnum '_fixtime_coeff_15.5.mat']), 'wbody');
    data2 = load(fullfile(rootdir, 'BodyPart_GazeTime_LabPhone', [vdnum '_fixtime_coeff_28.5.mat']),'wbody');
    
    combined = [data1.wbody; data2.wbody];
    bodyparts_collect.(vdname) = mean(combined,1);
end


%% LabTobii gaze to whole body and non-body (background)

results_dir = fullfile(rootdir, 'BodyPart_GazeTime_Tobii');

summary_table = {};

for vd=1:numel(videos)
    vdnum = videos{vd};
    data = load(fullfile(results_dir, [vdnum '_fixtime_coeff_15.5.mat']));
    
    subjects = cellstr(data.subj_info);
    
    for s=1:numel(subjects)
        subj = subjects{s};

        if ismember(subj, idmatch_v2.EarlyID) 
            idx = find(strcmp(idmatch_v2.EarlyID, subj));
            subj_matched = idmatch_v2.RealID{idx};
            subjects{s} = subj_matched;
        else
            subj_matched = subj;
        end
    end
    
    grouplabels = {};
    for g=1:length(subjects)
        if data.subj_groups(g) == 1
            grouplabels{g} = 'A';
        else
            grouplabels{g} = 'C';
        end
    end
    
    vdname = ['mv' vdnum];
    frames_w_bodyparts = bodyparts_collect.(vdname); 
    % if mean is 0, it most likely means this frame didn't have any body features
    
    onscreen = mean(data.onscreen,2);
    for g=1:length(subjects)
        subj = subjects{g};
        group = grouplabels{g};
        
        temp = cell(1,5);
        temp{1} = subj;
        temp{2} = group;
        temp{3} = vdnum;
        temp{4} = 'onscreen';
        temp{5} = onscreen(g);
        summary_table = [summary_table; temp];
    end

    for ft=1:numel(bodyfeatures)
        feat = bodyfeatures{ft};
        sumNorm = [];
        
        data_onscreen = int8(data.(feat));
        data_onscreen(~data.onscreen) = -1;  % code frames with off-screen gaze as -1       
        data_selected = data_onscreen(:, frames_w_bodyparts~=0); % selected frames with body features only

        sumNorm = sum(data_selected == 1,2) ./ sum(data_selected ~= -1,2);
        
        for g=1:length(subjects)
            subj = subjects{g};
            group = grouplabels{g};
            lookingtime = sumNorm(g);
            
            temp = cell(1,5);
            temp{1} = subj;
            temp{2} = group;
            temp{3} = vdnum;
            temp{4} = feat;
            temp{5} = lookingtime;
            
            summary_table = [summary_table; temp];
        end
    end
    
    clear data
end

T = cell2table(summary_table, "VariableNames", ["subj_id" "grouplabel" "video" "feature" "lookingtime" ]);
writetable(T, fullfile(results_dir, 'labtobii_bodyparts_summary.csv')) 


%% LabTobii gaze to face

results_dir = fullfile(rootdir, 'Face_GazeTime_Tobii');

summary_table = {};

for vd=1:numel(videos)
    vdnum = videos{vd};
    data = load(fullfile(results_dir, [vdnum '_fixtime_coeff_15.5.mat']));
    
    subjects = cellstr(data.subjs);
    
    for s=1:numel(subjects)
        subj = subjects{s};

        if ismember(subj, idmatch_v2.EarlyID) 
            idx = find(strcmp(idmatch_v2.EarlyID, subj));
            subj_matched = idmatch_v2.RealID{idx};
            subjects{s} = subj_matched;
        else
            subj_matched = subj;
        end
    end
    
    grouplabels = {};
    for g=1:length(subjects)
        if data.subj_groups(g) == 1
            grouplabels{g} = 'A';
        else
            grouplabels{g} = 'C';
        end
    end
    
    sumNorm = [];

    for ft=1:numel(facefeatures)
        feat = facefeatures{ft};
        
        % saving non-normalized data is not necessary here (-1 includes frames that have no face)
        sumNorm.(feat) = sum(data.(feat) == 1,2) ./ sum(data.(feat) ~= -1,2);
        
        for g=1:length(subjects)
            subj = subjects{g};
            group = grouplabels{g};
            lookingtime = sumNorm.(feat)(g);
            
            temp = cell(1,5);
            temp{1} = subj;
            temp{2} = group;
            temp{3} = vdnum;
            temp{4} = feat;
            temp{5} = lookingtime;
            
            summary_table = [summary_table; temp];
        end
    end
    
    clear data
end

T = cell2table(summary_table, "VariableNames", ["subj_id" "grouplabel" "video" "feature" "lookingtime"]);
writetable(T, fullfile(results_dir, 'labtobii_face_summary.csv')) 


%% LabPhone gaze to whole body and non-body (background)

results_dir = fullfile(rootdir, 'BodyPart_GazeTime_LabPhone');

summary_table = {};

for vd=1:numel(videos)
    vdnum = videos{vd};
    data = load(fullfile(results_dir, [vdnum '_fixtime_coeff_28.5.mat']));
    
    subjects = cellstr(data.subj_info);
    grouplabels = {};
    for g=1:length(subjects)
        if data.subj_groups(g) == 1
            grouplabels{g} = 'A';
        else
            grouplabels{g} = 'C';
        end
    end
    
    vdname = ['mv' vdnum];
    frames_w_bodyparts = bodyparts_collect.(vdname); 
    % if mean is 0, it most likely means this frame didn't have any body features
    
    onscreen = mean(data.onscreen,2);

    for ft=1:numel(bodyfeatures)
        feat = bodyfeatures{ft};
        sumNorm = [];
        
        data_onscreen = int8(data.(feat));
        data_onscreen(~data.onscreen) = -1;  % code frames with off-screen gaze as -1       
        data_selected = data_onscreen(:, frames_w_bodyparts~=0); % selected frames with body features only

        sumNorm = sum(data_selected == 1,2) ./ sum(data_selected ~= -1,2);
        
        for g=1:length(subjects)
            subj = subjects{g};
            group = grouplabels{g};
            lookingtime = sumNorm(g);
            
            temp = cell(1,5);
            temp{1} = subj;
            temp{2} = group;
            temp{3} = vdnum;
            temp{4} = feat;
            temp{5} = lookingtime;
            
            summary_table = [summary_table; temp];
        end
    end
    
    clear data
end

T = cell2table(summary_table, "VariableNames", ["subj_id" "grouplabel" "video" "feature" "lookingtime" ]);
writetable(T, fullfile(results_dir, 'labphone_bodyparts_summary.csv')) 


%% LabPhone gaze to face

results_dir = fullfile(rootdir, 'Face_GazeTime_LabPhone');

summary_table = {};

for vd=1:numel(videos)
    vdnum = videos{vd};
    data = load(fullfile(results_dir, [vdnum '_fixtime_coeff_28.5.mat']));
    
    subjects = cellstr(data.subjs);
    
    for s=1:numel(subjects)
        subj = subjects{s};

        if ismember(subj, idmatch_v2.EarlyID) 
            idx = find(strcmp(idmatch_v2.EarlyID, subj));
            subj_matched = idmatch_v2.RealID{idx};
            subjects{s} = subj_matched;
        else
            subj_matched = subj;
        end
    end
    
    grouplabels = {};
    for g=1:length(subjects)
        if data.subj_groups(g) == 1
            grouplabels{g} = 'A';
        else
            grouplabels{g} = 'C';
        end
    end
    

    sumNorm = [];

    for ft=1:numel(facefeatures)
        feat = facefeatures{ft};

        sumNorm.(feat) = sum(data.(feat) == 1,2) ./ sum(data.(feat) ~= -1,2);
        
        for g=1:length(subjects)
            subj = subjects{g};
            group = grouplabels{g};
            lookingtime = sumNorm.(feat)(g);
            
            temp = cell(1,5);
            temp{1} = subj;
            temp{2} = group;
            temp{3} = vdnum;
            temp{4} = feat;
            temp{5} = lookingtime;
            
            summary_table = [summary_table; temp];
        end
    end
    
    clear data
end

T = cell2table(summary_table, "VariableNames", ["subj_id" "grouplabel" "video" "feature" "lookingtime" ]);
writetable(T, fullfile(results_dir, 'labphone_face_summary.csv')) 

