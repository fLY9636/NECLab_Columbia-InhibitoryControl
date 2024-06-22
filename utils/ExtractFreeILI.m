function subj_ILI = ExtractFreeILI(loopIDX, MetaData, OFCIDX, ANIMAL_IDs, ANIMAL_VARs, MetaData_files, freewin_sz)
tic
subj_ILI = cell(1, length(loopIDX));
% freewin_sz = 4;
% ANIMAL_IDs = ANIMAL_IDs(loopIDX);
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    IDX = ANIMAL_VARs.(ANIMAL_IDs{i});
    session_ILI = [];
    for j = IDX(:).'
        m = load(MetaData_files{j});
        for k = 2:height(m.MetaData)
            if isempty(m.MetaData.Individual_Licks{k}) == 0 && isempty(m.MetaData.Individual_Licks{k-1}) == 0
                idx1 = [m.MetaData.Reward_Onset(k-1) m.MetaData.Punish_Onset(k-1)];
                boundary1 = idx1(~isnan(idx1))+freewin_sz;
                boundary2 = m.MetaData.Tone_Onset(k);
                % boundary-included ILI
                gap_lickraster = [boundary1; m.MetaData.Individual_Licks{k-1}(find(m.MetaData.Individual_Licks{k-1}(:, 1)>=boundary1), 1); ...
                    m.MetaData.Individual_Licks{k}(find(m.MetaData.Individual_Licks{k}(:, 1)<=boundary2), 1); boundary2];
%                 gap_lickraster = [boundary1; m.MetaData.Individual_Licks{k-1}(find(m.MetaData.Individual_Licks{k-1}(:, 1)>=boundary1), 1); ...
%                         m.MetaData.Individual_Licks{k}(find(m.MetaData.Individual_Licks{k}(:, 1)<=boundary2), 1)];
                % boundary-less ILI
%                 gap_lickraster = [m.MetaData.Individual_Licks{k-1}(find(m.MetaData.Individual_Licks{k-1}(:, 1)>=boundary1), 1); ...
%                     m.MetaData.Individual_Licks{k}(find(m.MetaData.Individual_Licks{k}(:, 1)<=boundary2), 1)];
                gap_lickraster = gap_lickraster-m.MetaData.Tone_Onset(k);
                if length(gap_lickraster)>2
                    session_ILI = [session_ILI; diff(gap_lickraster(1:end-1))];
                end
                % semi-boundary no-lick-period cutoff X[.....|]
%                 if length(gap_lickraster)==2
%                     session_ILI = [session_ILI; diff(gap_lickraster)];
%                 end
                % semi-boundary 
%                 if length(gap_lickraster)==1
%                     session_ILI(end+1, 1) = (boundary2-boundary1);
%                 end
%                 if length(gap_lickraster)==2
%                     session_ILI = [session_ILI; diff(gap_lickraster)];
%                 end
                % semi-boundary one-lick-period cutoff X[...|...|]
%                 if length(gap_lickraster)==1
%                     session_ILI(end+1, 1) = (boundary2-boundary1)/2;
%                 end
                % semi-boundary one-lick period halved 
                if length(gap_lickraster)==2
                    session_ILI(end+1, 1) = (boundary2-boundary1);
                end
                if length(gap_lickraster)==3
                    session_ILI = [session_ILI; min([diff(gap_lickraster(1:end-1)) diff(gap_lickraster(2:end))])];
                end
                if length(gap_lickraster)>3
                    session_ILI = [session_ILI; diff(gap_lickraster(2:end-1))];
                end
%                     session_ILI(end+1, 1) = (boundary2-boundary1)/2;
            end
        end
    end
    subj_ILI{cnt} = session_ILI;
end
return
toc