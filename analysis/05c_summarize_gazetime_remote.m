
rootdir = 'results';

videos.session1  = {'LXX01'; 'L0102'; 'L0103'; 'LXX04'; 'L0105'; 'L0106'; 'L0107'; 'L0108'; 'L0109'};
videos.session2  = {'LXX01'; 'L0202'; 'L0203'; 'LXX04'; 'L0205'; 'L0206'; 'L0207'; 'L0208'; 'L0209'};
videos.session3  = {'LXX01'; 'L0302'; 'L0303'; 'LXX04'; 'L0305'; 'L0306'; 'L0307'; 'L0308'; 'L0309'};
videos.session4  = {'LXX01'; 'L0402'; 'L0403'; 'LXX04'; 'L0405'; 'L0406'; 'L0407'; 'L0408'; 'L0409'};
videos.session5  = {'LXX01'; 'L0502'; 'L0503'; 'LXX04'; 'L0505'; 'L0506'; 'L0507'; 'L0508'; 'L0509'};
videos.session6  = {'LXX01'; 'L0602'; 'L0603'; 'LXX04'; 'L0605'; 'L0606'; 'L0607'; 'L0608'; 'L0609'};
videos.session7  = {'LXX01'; 'L0702'; 'L0703'; 'LXX04'; 'L0705'; 'L0706'; 'L0707'; 'L0708'; 'L0709'};
videos.session8  = {'LXX01'; 'L0802'; 'L0803'; 'LXX04'; 'L0805'; 'L0806'; 'L0807'; 'L0808'; 'L0809'};
videos.session9  = {'LXX01'; 'L0902'; 'L0903'; 'LXX04'; 'L0905'; 'L0906'; 'L0907'; 'L0908'; 'L0909'};
videos.session10 = {'LXX01'; 'L1002'; 'L1003'; 'LXX04'; 'L1005'; 'L1006'; 'L1007'; 'L1008'; 'L1009'};

bodyfeatures = {'wbody'; 'nonbody'};
facefeatures = {'face_gaze'; 'center_gaze'};


%% index video frames without any human body parts 

% assumption: when no one has gaze to wbody features across all participants, it means there's no body feature on the frame. 

bodyparts_collect = [];
for ses=1:10
    sessionName = ['session' num2str(ses)];
    this_videos = videos.(sessionName); 

    for vd=1:numel(this_videos)
        vdnum = this_videos{vd};
        vdname = ['mv' vdnum];
        data2 = load(fullfile(rootdir, ['BodyPart_GazeTime_RemotePhone_ses' num2str(ses)], [vdnum '_fixtime_coeff_28.5.mat']),'wbody');

        combined = data2.wbody;
        bodyparts_collect.(sessionName).(vdname) = mean(combined,1);

    end
end

%% RemotePhone gaze to whole body and non-body (background)

for ses=1:10
    sessionName = ['session' num2str(ses)];
    this_videos = videos.(sessionName); 
    
    results_dir = fullfile(rootdir, ['BodyPart_GazeTime_RemotePhone_ses' num2str(ses)]);

    summary_table = {};
    
    for vd=1:numel(this_videos)
        vdnum = this_videos{vd};
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
        frames_w_bodyparts = bodyparts_collect.(sessionName).(vdname); 
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

                temp = cell(1,6);
                temp{1} = subj;
                temp{2} = group;
                temp{3} = sessionName;
                temp{4} = vdnum;
                temp{5} = feat;
                temp{6} = lookingtime;

                summary_table = [summary_table; temp];
            end
        end

        clear data
    end

    T = cell2table(summary_table, "VariableNames", ["subj_id" "grouplabel" "session" "video" "feature" "lookingtime" ]);
    writetable(T, fullfile(results_dir, ['remotephone_' sessionName '_bodyparts_summary.csv']))     
    
end






%% RemotePhone gaze to face


for ses=1:10
    sessionName = ['session' num2str(ses)];
    this_videos = videos.(sessionName); 
    
    results_dir = fullfile(rootdir, ['Face_GazeTime_RemotePhone_ses' num2str(ses)]);

    summary_table = {};
    
    for vd=1:numel(this_videos)
        vdnum = this_videos{vd};
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

                temp = cell(1,6);
                temp{1} = subj;
                temp{2} = group;
                temp{3} = sessionName;
                temp{4} = vdnum;
                temp{5} = feat;
                temp{6} = lookingtime;

                summary_table = [summary_table; temp];
            end
        end

        clear data
    end

    T = cell2table(summary_table, "VariableNames", ["subj_id" "grouplabel" "session" "video" "feature" "lookingtime" ]);
    writetable(T, fullfile(results_dir, ['remotephone_' sessionName '_face_summary.csv']))     
    
    
end
 