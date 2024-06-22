[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 0, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);

num_bins = 10;
session_AUC_NE = NaN(len, 1);
session_AUC_Ach = NaN(len, 1);
session_AUC_NA = NaN(len, 1);
session_ROC_NE = cell(len, 1);
session_ROC_Ach = cell(len, 1);
session_ROC_NA = cell(len, 1);
WINst = 2*120+1;
WINed = 4*120;

% for i = 1:length(OFCIDX)
%     dF_NE = nanmean(session_NE_atOutcome_P{i}(:, WINst:WINed), 2);
%     dS_NE = nanmean(session_NE_atOutcome_R{i}(:, WINst:WINed), 2);
%     dF_Ach = nanmean(session_Ach_atOutcome_P{i}(:, WINst:WINed), 2);
%     dS_Ach = nanmean(session_Ach_atOutcome_R{i}(:, WINst:WINed), 2);
%     dF_NA = nanmean(session_NA_atOutcome_P{i}(:, WINst:WINed), 2);
%     dS_NA = nanmean(session_NA_atOutcome_R{i}(:, WINst:WINed), 2);
%     thresholds_NE = min([dF_NE; dS_NE]):(max([dF_NE; dS_NE])-min([dF_NE; dS_NE]))/num_bins:max([dF_NE; dS_NE]);
%     [p2, AUC] = CalcROC(dF_NE, dS_NE, thresholds_NE);
%     session_AUC_NE(i) = abs(AUC-0.5);
%     if AUC<0.5
%         session_ROC_NE{i} = p2([2 1], :);
%     else
%         session_ROC_NE{i} = p2;
%     end
%     thresholds_Ach = min([dF_Ach; dS_Ach]):(max([dF_Ach; dS_Ach])-min([dF_Ach; dS_Ach]))/num_bins:max([dF_Ach; dS_Ach]);
%     [p2, AUC] = CalcROC(dF_Ach, dS_Ach, thresholds_Ach);
%     session_AUC_Ach(i) = abs(AUC-0.5);
%     if AUC<0.5
%         session_ROC_Ach{i} = p2([2 1], :);
%     else
%         session_ROC_Ach{i} = p2;
%     end
%     thresholds_NA = min([dF_NA; dS_NA]):(max([dF_NA; dS_NA])-min([dF_NA; dS_NA]))/num_bins:max([dF_NA; dS_NA]);
%     [p2, AUC] = CalcROC(dF_NA, dS_NA, thresholds_NA);
%     session_AUC_NA(i) = abs(AUC-0.5);
%     if AUC<0.5
%         session_ROC_NA{i} = p2([2 1], :);
%     else
%         session_ROC_NA{i} = p2;
%     end
% end

num_wins = 4;
winsz = (WINed-WINst)/num_wins;
for i = 1:length(OFCIDX)
    win_AUC_NE = NaN(1, num_wins);
    win_AUC_Ach = NaN(1, num_wins); 
    win_AUC_NA = NaN(1, num_wins); 
    for j = 1:num_wins
        winst = WINst+winsz*(j-1); wined = winst+winsz-1;
        dF_NE = nanmean(session_NE_atOutcome_P{i}(:, winst:wined), 2);
        dS_NE = nanmean(session_NE_atOutcome_R{i}(:, winst:wined), 2);
        dF_Ach = nanmean(session_Ach_atOutcome_P{i}(:, winst:wined), 2);
        dS_Ach = nanmean(session_Ach_atOutcome_R{i}(:, winst:wined), 2);
        dF_NA = session_sync_LDA_P{i};
        dS_NA = session_sync_LDA_R{i};
        thresholds_NE = min([dF_NE; dS_NE]):(max([dF_NE; dS_NE])-min([dF_NE; dS_NE]))/num_bins:max([dF_NE; dS_NE]);
        [p2, AUC] = CalcROC(dF_NE, dS_NE, thresholds_NE);
        win_AUC_NE(j) = abs(AUC-0.5);
        thresholds_Ach = min([dF_Ach; dS_Ach]):(max([dF_Ach; dS_Ach])-min([dF_Ach; dS_Ach]))/num_bins:max([dF_Ach; dS_Ach]);
        [p2, AUC] = CalcROC(dF_Ach, dS_Ach, thresholds_Ach);
        win_AUC_Ach(j) = abs(AUC-0.5);
        thresholds_NA = min([dF_NA; dS_NA]):(max([dF_NA; dS_NA])-min([dF_NA; dS_NA]))/num_bins:max([dF_NA; dS_NA]);
        [p2, AUC] = CalcROC(dF_NA, dS_NA, thresholds_NA);
        win_AUC_NA(j) = abs(AUC-0.5);
    end
    session_AUC_NE(i) = mean(win_AUC_NE);
    session_AUC_Ach(i) = mean(win_AUC_Ach);
    session_AUC_NA(i) = mean(win_AUC_NA);
end
%% Saline AUC Lumped
subjMEAN_AUC_NE = [];
subjMEAN_AUC_Ach = [];
subjMEAN_AUC_NA = [];
for i = 1:length(subjIDX)
    subjMEAN_AUC_NE(end+1, 1) = mean(session_AUC_NE(subjIDX{i}));
    subjMEAN_AUC_Ach(end+1, 1) = mean(session_AUC_Ach(subjIDX{i}));
    subjMEAN_AUC_NA(end+1, 1) = mean(session_AUC_NA(subjIDX{i}));
end
%%
figure; hold on
bar(1, mean(subjMEAN_AUC_NE))
bar(2, mean(subjMEAN_AUC_Ach))
bar(3, mean(subjMEAN_AUC_NA))
for i = 1:length(subjIDX)
    plot([1 2 3], [subjMEAN_AUC_NE(i) subjMEAN_AUC_Ach(i) subjMEAN_AUC_NA(i)], '-ok')
end
errorbar([1 2 3], [mean(subjMEAN_AUC_NE) mean(subjMEAN_AUC_Ach) mean(subjMEAN_AUC_Ach)], ...
    [SEM(subjMEAN_AUC_NE) SEM(subjMEAN_AUC_Ach) SEM(subjMEAN_AUC_Ach)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('|AUROC-0.5| (a.u.)')

%% Saline AUC Lumped
subjMEAN_Saline_AUC_NE = [];
subjMEAN_Saline_AUC_Ach = [];
subjMEAN_Saline_AUC_NA = [];
for i = 1:length(subjIDX)
    subjMEAN_Saline_AUC_NE(end+1, 1) = mean(session_AUC_NE(subjIDX{i}));
    subjMEAN_Saline_AUC_Ach(end+1, 1) = mean(session_AUC_Ach(subjIDX{i}));
    subjMEAN_Saline_AUC_NA(end+1, 1) = mean(session_AUC_NA(subjIDX{i}));
end
sessMEAN_Saline_AUC_NE = session_AUC_NE;
sessMEAN_Saline_AUC_Ach = session_AUC_Ach;
sessMEAN_Saline_AUC_NA = session_AUC_NA;
%% CNO AUC Lumped
subjMEAN_CNO_AUC_NE = [];
subjMEAN_CNO_AUC_Ach = [];
subjMEAN_CNO_AUC_NA = [];
for i = 1:length(subjIDX)
    subjMEAN_CNO_AUC_NE(end+1, 1) = mean(session_AUC_NE(subjIDX{i}));
    subjMEAN_CNO_AUC_Ach(end+1, 1) = mean(session_AUC_Ach(subjIDX{i}));
    subjMEAN_CNO_AUC_NA(end+1, 1) = mean(session_AUC_NA(subjIDX{i}));
end
sessMEAN_CNO_AUC_NE = session_AUC_NE;
sessMEAN_CNO_AUC_Ach = session_AUC_Ach;
sessMEAN_CNO_AUC_NA = session_AUC_NA;
%%
figure; hold on
bar(1, mean(subjMEAN_Saline_AUC_NE))
bar(2, mean(subjMEAN_CNO_AUC_NE))
bar(4, mean(subjMEAN_Saline_AUC_Ach))
bar(5, mean(subjMEAN_CNO_AUC_Ach))
bar(7, mean(subjMEAN_Saline_AUC_NA))
bar(8, mean(subjMEAN_CNO_AUC_NA))
for i = 1:length(subjMEAN_Saline_AUC_NE)
    plot(1:2, [subjMEAN_Saline_AUC_NE(i) subjMEAN_CNO_AUC_NE(i)], '-ok')
    plot(4:5, [subjMEAN_Saline_AUC_Ach(i) subjMEAN_CNO_AUC_Ach(i)], '-ok')
    plot(7:8, [subjMEAN_Saline_AUC_NA(i) subjMEAN_CNO_AUC_NA(i)], '-ok')
end
errorbar([1 2 4 5 7 8], [mean(subjMEAN_Saline_AUC_NE) mean(subjMEAN_CNO_AUC_NE) mean(subjMEAN_Saline_AUC_Ach) mean(subjMEAN_CNO_AUC_Ach) ...
    mean(subjMEAN_Saline_AUC_NA) mean(subjMEAN_CNO_AUC_NA)],...
     [SEM(subjMEAN_Saline_AUC_NE) SEM(subjMEAN_CNO_AUC_NE) SEM(subjMEAN_Saline_AUC_Ach) SEM(subjMEAN_CNO_AUC_Ach) ...
     SEM(subjMEAN_Saline_AUC_NA) SEM(subjMEAN_CNO_AUC_NA)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('|AUROC-0.5| (a.u.)')

% figure; hold on
% bar(1, mean(subjMEAN_Saline_AUC_NA))
% bar(2, mean(subjMEAN_CNO_AUC_NA))
% for i = 1:length(subjMEAN_Saline_AUC_NA)
%     plot(1:2, [subjMEAN_Saline_AUC_NA(i) subjMEAN_CNO_AUC_NA(i)], '-ok')
% end
% errorbar([1 2], [mean(subjMEAN_Saline_AUC_NA) mean(subjMEAN_CNO_AUC_NA)], [SEM(subjMEAN_Saline_AUC_NA) SEM(subjMEAN_CNO_AUC_NA)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('\Delta AUROC (a.u.)')
%%
figure; hold on
bar(1, mean(sessMEAN_Saline_AUC_NE))
bar(2, mean(sessMEAN_CNO_AUC_NE))
bar(4, mean(sessMEAN_Saline_AUC_Ach))
bar(5, mean(sessMEAN_CNO_AUC_Ach))
scatter(ones(length(sessMEAN_Saline_AUC_NE), 1), sessMEAN_Saline_AUC_NE)
scatter(2*ones(length(sessMEAN_CNO_AUC_NE), 1), sessMEAN_CNO_AUC_NE)
scatter(4*ones(length(sessMEAN_Saline_AUC_Ach), 1), sessMEAN_Saline_AUC_Ach)
scatter(5*ones(length(sessMEAN_CNO_AUC_Ach), 1), sessMEAN_CNO_AUC_Ach)
errorbar([1 2 4 5], [mean(sessMEAN_Saline_AUC_NE) mean(sessMEAN_CNO_AUC_NE) mean(sessMEAN_Saline_AUC_Ach) mean(sessMEAN_CNO_AUC_Ach)],...
     [SEM(sessMEAN_Saline_AUC_NE) SEM(sessMEAN_CNO_AUC_NE) SEM(sessMEAN_Saline_AUC_Ach) SEM(sessMEAN_CNO_AUC_Ach)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('\Delta AUROC (a.u.)')

% figure; hold on
% bar(1, mean(sessMEAN_Saline_AUC_NA))
% bar(2, mean(sessMEAN_CNO_AUC_NA))
% scatter(ones(length(sessMEAN_Saline_AUC_NE), 1), sessMEAN_Saline_AUC_NE)
% scatter(2*ones(length(sessMEAN_CNO_AUC_NE), 1), sessMEAN_CNO_AUC_NE)
% errorbar([1 2], [mean(sessMEAN_Saline_AUC_NA) mean(sessMEAN_CNO_AUC_NA)], [SEM(sessMEAN_CNO_AUC_NE) SEM(sessMEAN_CNO_AUC_NE)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('|AUROC-0,5| (a.u.)')



