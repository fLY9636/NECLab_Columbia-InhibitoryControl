%% Identify session number of treatment groups
% CLONIDINE_DOSEs = {'C', 'L', 'H'};
% CLONIDINE_VARs = struct; 
% for i = 1:length(CLONIDINE_DOSEs)
%     CLONIDINE_VARs.(CLONIDINE_DOSEs{i}) = [];
% end
% for i = 1:length(Behavior_files)
%     b = load(Behavior_files{i});
%     if strcmp(b.Data.Memo, 'C')==1
%         CLONIDINE_VARs.C(end+1) = i;
%     end
%     if strcmp(b.Data.Memo, 'L')==1
%         CLONIDINE_VARs.L(end+1) = i;
%     end
%     if strcmp(b.Data.Memo, 'H')==1
%         CLONIDINE_VARs.H(end+1) = i;
%     end 
% end  
% %% Calculate and visulize withholding success rate of each treatment group
% Lumped_treatment_subj_performance = {};
% for i = 1:13
%     Lumped_treatment_subj_performance = [Lumped_treatment_subj_performance; {[] [] []}];
% end
% RT_threshold = 1.5; % set reaction time threshold
% for i = 1:length(ANIMAL_IDs)
%     for j = ANIMAL_VARs.(ANIMAL_IDs{i})
%         m = load(MetaData_files{j});
%         b = load(Behavior_files{j});
%         num_Rwd = 0;
%         for k = 1:height(m.MetaData)
%             if isempty(m.MetaData.Individual_Licks{k})==0
%                 if isempty(find(m.MetaData.Individual_Licks{k}(:, 1)>=m.MetaData.Reward_Onset(k)&m.MetaData.Individual_Licks{k}(:, 1)<=m.MetaData.Reward_Onset(k)+RT_threshold))==0
%                     num_Rwd = num_Rwd+1;
%                 end
%             end
%         end
%         if strcmp(b.Data.Memo, 'L')==1
%             Lumped_treatment_subj_performance{i, 2}(end+1, 1) = num_Rwd/height(m.MetaData);
%         end
%         if strcmp(b.Data.Memo, 'H')==1
%             Lumped_treatment_subj_performance{i, 3}(end+1, 1) = num_Rwd/height(m.MetaData);
%         end
%         if strcmp(b.Data.Memo, 'C')==1
%             Lumped_treatment_subj_performance{i, 1}(end+1, 1) = num_Rwd/height(m.MetaData);
%         end
%     end
% end
%%
[Behavior_files, Phot_files, Pupil_files, MetaData_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc(ROOTDIR, 6, 0);
%% Calculate and visulize withholding success rate of non-treatment group
loopID = [3:5]; 
% sessionID = {2:4 3:5};
session_noCNO = find(MetaData_DRREADs1.CNO==0); session_CNO = find(MetaData_DRREADs1.CNO==1);
% session_LoCNO = [5:7 9:10 25:29]'; session_HiCNO = setdiff(session_CNO, session_LoCNO);
% session_CNO = find(MetaData.CNO==1); session_noCNO = find(MetaData.CNO==0); 
Lumped_subj_performance_duration = cell(1, length(loopID));
Lumped_subj_performance = cell(1, length(loopID));
Lumped_subj_punish = cell(1, length(loopID));
Lumped_subj_slowresp = cell(1, length(loopID));
Lumped_subj_RT = cell(1, length(loopID));
Lumped_subj_sessionRT = cell(1, length(loopID));
RT_threshold = 1.5; % set reaction time threshold
cnt = 0;
for i = loopID
    cnt = cnt+1;
%     for j = ANIMAL_VARs.(ANIMAL_IDs{i})
%         if ismember(j, ([1:5 9:13 17:21]))==0
%             continue;
%         end
%     for j = (intersect(ANIMAL_VARs.(ANIMAL_IDs{i}), session_noCNO))'
%         if ~ismember(j, ([46, 48, 49, 51, 54, 60, 61, 63, 65, 66, 75, 77, 79, 80, 81]))
%         if ~ismember(j, ([3, 5, 7, 9, 10, 18, 19, 21, 24, 25]))
%             continue;
%         end
    for j = (intersect(ANIMAL_VARs.(ANIMAL_IDs{i}), session_CNO))'
        m = load(MetaData_files{j});
        b = load(Behavior_files{j});
        num_Rwd = 0; num_punish = 0; num_slowresp = 0;
        num_Rwd_short = 0; num_Rwd_medium = 0; num_Rwd_long = 0;
        num_short = 0; num_medium = 0; num_long = 0;
        sessionRT = [];
        for k = 1:height(m.MetaData)
            if ~isempty(m.MetaData.Individual_Licks{k})
                if isnan(m.MetaData.Punish_Onset(k))
                    if ~isempty(find(m.MetaData.Individual_Licks{k}(:, 1)>=m.MetaData.Reward_Onset(k)&m.MetaData.Individual_Licks{k}(:, 1)<=m.MetaData.Reward_Onset(k)+RT_threshold))
                        num_Rwd = num_Rwd+1;
                    end
                    if isempty(find(m.MetaData.Individual_Licks{k}(:, 1)>=m.MetaData.Reward_Onset(k)&m.MetaData.Individual_Licks{k}(:, 1)<=m.MetaData.Reward_Onset(k)+RT_threshold))
                        num_slowresp = num_slowresp+1;
                    end
                    if ~isempty(find(m.MetaData.Individual_Licks{k}(:, 1)>=m.MetaData.Reward_Onset(k)))
                        reactlick = find(m.MetaData.Individual_Licks{k}(:, 1)>=m.MetaData.Reward_Onset(k));
                        sessionRT(end+1, 1) = m.MetaData.Individual_Licks{k}(reactlick(1), 1)-m.MetaData.Reward_Onset(k);
                        Lumped_subj_RT{cnt}(end+1, 1) = m.MetaData.Individual_Licks{k}(reactlick(1), 1)-m.MetaData.Reward_Onset(k);
                    else
                        sessionRT(end+1, 1) = -1;
                        Lumped_subj_RT{cnt}(end+1, 1) = -1;
                    end
                end
                if ~isnan(m.MetaData.Punish_Onset(k))
                    num_punish = num_punish+1;
                end
            end
            if isempty(m.MetaData.Individual_Licks{k})
                num_slowresp = num_slowresp+1;
                sessionRT(end+1, 1) = -1;
                Lumped_subj_RT{cnt}(end+1, 1) = -1;
            end
        end
        Lumped_subj_performance{cnt}(end+1, 1) = num_Rwd/height(m.MetaData);
        Lumped_subj_punish{cnt}(end+1, 1) = num_punish/height(m.MetaData);
        Lumped_subj_slowresp{cnt}(end+1, 1) = num_slowresp/height(m.MetaData);
        Lumped_subj_sessionRT{cnt}(end+1, 1) = nanmean(sessionRT(find(sessionRT>=0)));
    end
end
%%
Lumped_majority_subj_performance = Lumped_subj_performance;
Lumped_majority_subj_punish = Lumped_subj_punish;
% Lumped_majority_subj_slowresp = cellfun(@plus, Lumped_majority_subj_performance, Lumped_majority_subj_punish, 'UniformOutput', false);
Lumped_majority_subj_punish_durations = Lumped_subj_punish_durations;
Lumped_majority_subj_performance_bins = Lumped_subj_performance_bins;
Lumped_majority_subj_slowresp = Lumped_subj_slowresp;
Lumped_majority_subj_RT = Lumped_subj_RT;
Lumped_majority_subj_sessionRT = Lumped_subj_sessionRT;
Lumped_majority_subj_performance_duration = Lumped_subj_performance_duration;
Lumped_majority_subj_failATduration = Lumped_subj_failATduration;
% Lumped_majority_subj_impulsive_ratio = Lumped_subj_impulsive_ratio;

%%
Lumped_treatment_subj_performance = Lumped_subj_performance;
Lumped_treatment_subj_punish = Lumped_subj_punish;
% Lumped_treatment_subj_slowresp = cellfun(@plus, Lumped_treatment_subj_performance, Lumped_treatment_subj_punish, 'UniformOutput', false);
Lumped_treatment_subj_punish_durations = Lumped_subj_punish_durations;
Lumped_treatment_subj_performance_bins = Lumped_subj_performance_bins;
Lumped_treatment_subj_slowresp = Lumped_subj_slowresp;
Lumped_treatment_subj_RT = Lumped_subj_RT;
Lumped_treatment_subj_sessionRT = Lumped_subj_sessionRT;
Lumped_treatment_subj_performance_duration = Lumped_subj_performance_duration;
Lumped_treatment_subj_failATduration = Lumped_subj_failATduration;
% Lumped_treatment_subj_impulsive_ratio = Lumped_subj_impulsive_ratio;
%%
figure
centers = 0.5:1:11.5;
ecdf(cell2mat(Lumped_majority_subj_punish_durations'))
hold on
ecdf(cell2mat(Lumped_treatment_subj_punish_durations'))
%%
num_duration_bins = 7;
duration_bins = 5:(12-5)/num_duration_bins:12;
bin_config = {[] [] [] [] [] [] []};
durationbin_failprogress_majority = bin_config; durationbin_failprogress_treatment = bin_config;
majority_set = cell2mat(Lumped_majority_subj_failATduration'); treatment_set = cell2mat(Lumped_treatment_subj_failATduration');
for iter = 1:num_duration_bins
    durationbin_failprogress_majority{iter} = [durationbin_failprogress_majority{iter}; majority_set(find(majority_set(:, 2)>duration_bins(iter)&majority_set(:, 2)<duration_bins(iter+1)), 1)];
    durationbin_failprogress_treatment{iter} = [durationbin_failprogress_treatment{iter}; treatment_set(find(treatment_set(:, 2)>duration_bins(iter)&treatment_set(:, 2)<duration_bins(iter+1)), 1)];
end
figure; hold on
MEAN1 = cellfun(@nanmean, durationbin_failprogress_majority); SEM1 = cellfun(@SEM, durationbin_failprogress_majority);
MEAN2 = cellfun(@nanmean, durationbin_failprogress_treatment); SEM2 = cellfun(@SEM, durationbin_failprogress_treatment);
MEAN = [MEAN1' MEAN2']; sem = [SEM1' SEM2'];
nbars = 2;
groupwidth = min(0.8, nbars/(nbars+1.5));
bar(MEAN)
for i = 1:nbars
    x = (1:num_duration_bins) - groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
    errorbar(x, MEAN(:, i), sem(:, i)/1, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
end
%%
input1 = cell2mat(Lumped_majority_subj_performance_bins'); input2 = cell2mat(Lumped_treatment_subj_performance_bins');
MEAN1 = nanmean(input1, 1); SEM1 = SEM(input1);
MEAN2 = nanmean(input2, 1); SEM2 = SEM(input2);
% input3 = session_AUC_5HT; MEAN3 = mean(session_AUC_5HT, 1); SEM3 = SEM(session_AUC_5HT);
MEAN = [MEAN1' MEAN2']; sem = [SEM1' SEM2'];
figure; hold on
nbars = 2;
groupwidth = min(0.8, nbars/(nbars+1.5));
bar(MEAN*100)
for i = 1:nbars
    x = (1:num_duration_bins) - groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
    errorbar(x, MEAN(:, i)*100, sem(:, i)/1*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
end
title('OFC')
ylabel('P( Outcome=Success | Tone Duration=t )   %')
xticks(1:7); xticklabels({'5~6', '6~7', '7~8', '8~9', '9~10', '10~11', '11~12'}); xtickangle(45)
xlabel('Tone Duration (s)')
title(['WoO = ' num2str(oppwinsz)])

num_duration_bins = 7;
duration_bins = 5:(12-5)/num_duration_bins:12;
bin_config = {[] [] [] [] [] [] []};
durationbin_RT_majority = bin_config; durationbin_RT_treatment = bin_config;
durationbin_Miss_majority = bin_config; durationbin_Miss_treatment = bin_config;
majority_set = cell2mat(Lumped_majority_subj_RT'); treatment_set = cell2mat(Lumped_treatment_RT');
for iter = 1:num_duration_bins
    IDX1 = find(majority_set(:, 4)>duration_bins(iter)&majority_set(:, 4)<duration_bins(iter+1));
    majority_temp = majority_set(IDX, 3); majority_temp = majority_temp(~isnan(majority_temp));
    durationbin_RT_majority{iter} = [durationbin_RT_majority{iter}; majority_set(IDX, 3)];
    IDX2 = find(treatment_set(:, 4)>duration_bins(iter)&treatment_set(:, 4)<duration_bins(iter+1));
    treatment_temp = treatment_set(IDX, 3); treatment_temp = treatment_temp(~isnan(treatment_temp));
    durationbin_RT_treatment{iter} = [durationbin_RT_treatment{iter}; treatment_set(IDX2, 3)];
end
%%
edges_failduration = 0:1:12; 
% input1 = cellfun(@(x) histcounts(x, edges_failduration, 'Normalization', 'probability'), Lumped_majority_subj_punish_durations', 'UniformOutput', false);
% input2 = cellfun(@(x) histcounts(x, edges_failduration, 'Normalization', 'probability'), Lumped_treatment_subj_punish_durations', 'UniformOutput', false);
% Lumped_MEAN1 = nanmean(cell2mat(input1), 1); Lumped_MEAN2 = nanmean(cell2mat(input2), 1); 
% Lumped_SEM1 = SEM(cell2mat(input1)); Lumped_SEM2 = SEM(cell2mat(input2));
% figure; hold on
% edges_failduration = edges_failduration(2:end);
% ShadedPlot(edges_failduration, Lumped_MEAN1, [0 0 0], 1, Lumped_SEM1, [0.9 0.8 0.7])
% plot(edges_failduration, Lumped_MEAN1, 'k', 'LineWidth', 1)
% ShadedPlot(edges_failduration, Lumped_MEAN2, [0 0.6 1], 1, Lumped_SEM2, [0.9 0.8 0.7])
% plot(edges_failduration, Lumped_MEAN2, 'Color', [0 0.6 1], 'LineWidth', 1)

input1 = cell2mat(Lumped_majority_subj_punish_durations');
input2 = cell2mat(Lumped_treatment_subj_punish_durations');
N1 = histcounts(input1, edges_failduration, 'Normalization', 'probability');
N2 = histcounts(input2, edges_failduration, 'Normalization', 'probability');
figure; hold on
edges_failduration = edges_failduration(2:end);
plot(edges_failduration, N1)
plot(edges_failduration, N2)
xlabel('Failed Withholding Duration (s)')

%%
input1 = cell2mat(Lumped_majority_subj_impulsive_ratio'); MEAN1 = nanmean(input1, 1); SEM1 = SEM(input1);
input2 = cell2mat(Lumped_treatment_subj_impulsive_ratio'); MEAN2 = nanmean(input2, 1); SEM2 = SEM(input2);
MEAN = [MEAN1' MEAN2']; sem = [SEM1' SEM2'];
figure; hold on
bar(1:2, MEAN*100)
errorbar(1:2, MEAN*100, sem/1*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
title('OFC')
ylabel('Impulsive ratio %')
%% Barplot of behavioral performance (control group and treatment group)
IDX = 1:28; subj_COLOR; h_plot = [];
%%%%%% comment if there is not a third dose %%%%%%%%%%%%%%%%%%%
input1 = Lumped_majority_subj_performance(IDX); input2 = Lumped_treatment_subj_performance(IDX);
input1 = input1(~cellfun(@isempty, input1)); input2 = input2(~cellfun(@isempty, input2));
figure; hold on
subplot(1, 3, 1); hold on
% bar([1 2], [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input2'), 1)]*100);
bar([1 2], [nanmean(cellfun(@nanmean, input1)) nanmean(cellfun(@nanmean, input2))]*100);
errorbar([1 2], [nanmean(cellfun(@nanmean, input1)) nanmean(cellfun(@nanmean, input2))]*100, ...
    [SEM(cell2mat(input1'))/1.1 SEM(cell2mat(input2'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
for i = 1:length(input1)
    h_plot(end+1) = plot([1 2], [nanmean(input1{i}) nanmean(input2{i})]*100, '-o', 'Color',  subj_COLOR{i}, 'MarkerEdgeColor', subj_COLOR{i})
end
ylim([0 100]); ylabel(['Successful withholding % \newline (lick within ' num2str(oppwinsz) 's after reward was delivered)']);
legend(h_plot, ANIMAL_IDs)
title(['WoO = ' num2str(oppwinsz)])
xticks([1 2]); xticklabels({'Saline', 'CNO'})
input1 = Lumped_majority_subj_punish(IDX); input2 = Lumped_treatment_subj_punish(IDX);
input1 = input1(~cellfun(@isempty, input1)); input2 = input2(~cellfun(@isempty, input2));
subplot(1, 3, 2); hold on
% bar([1 2], [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input2'), 1)]*100);
bar([1 2], [nanmean(cellfun(@nanmean, input1)) nanmean(cellfun(@nanmean, input2))]*100);
errorbar([1 2], [nanmean(cellfun(@nanmean, input1)) nanmean(cellfun(@nanmean, input2))]*100, ...
    [SEM(cell2mat(input1'))/1.1 SEM(cell2mat(input2'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
for i = 1:length(input1)
    cnt = cnt+1;
    plot([1 2], [nanmean(input1{i}) nanmean(input2{i})]*100, '-o', 'Color',  subj_COLOR{i}, 'MarkerEdgeColor', subj_COLOR{i})
end
ylim([0 100]); ylabel('Failed withholding %');
title(['WoO = ' num2str(oppwinsz)])
xticks([1 2]); xticklabels({'Saline', 'CNO'})
input1 = Lumped_majority_subj_slowresp(IDX); input2 = Lumped_treatment_subj_slowresp(IDX);
input1 = input1(~cellfun(@isempty, input1)); input2 = input2(~cellfun(@isempty, input2));
subplot(1, 3, 3); hold on
% bar(1:2, [nanmean(1-cell2mat(Lumped_majority_subj_punish')-cell2mat(Lumped_majority_subj_performance'), 1) ...
%     nanmean(1-cell2mat(Lumped_treatment_subj_punish')-cell2mat(Lumped_treatment_subj_performance'), 1)]*100);
% errorbar(1:2, [nanmean(1-cell2mat(Lumped_majority_subj_punish')-cell2mat(Lumped_majority_subj_performance'), 1) ...
%     nanmean(1-cell2mat(Lumped_treatment_subj_punish')-cell2mat(Lumped_treatment_subj_performance'), 1)]*100, ...
%     [SEM(1-cell2mat(Lumped_majority_subj_punish')-cell2mat(Lumped_majority_subj_performance'))/1.1 ...
%     SEM(1-cell2mat(Lumped_treatment_subj_punish')-cell2mat(Lumped_treatment_subj_performance'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% for i = 1:length(input1)
%     cnt = cnt+1;
%     plot(1:2, [nanmean(1-Lumped_majority_subj_performance{i}-Lumped_majority_subj_punish{i}) ...
%         nanmean(1-Lumped_treatment_subj_performance{i}-Lumped_treatment_subj_punish{i})]*100, '-o', 'Color',  [1 0.844 0], 'MarkerEdgeColor', [1 0.844 0])
% end
% bar([1 2], [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input2'), 1)]*100);
bar([1 2], [nanmean(cellfun(@nanmean, input1)) nanmean(cellfun(@nanmean, input2))]*100);
errorbar([1 2], [nanmean(cellfun(@nanmean, input1)) nanmean(cellfun(@nanmean, input2))]*100, ...
    [SEM(cell2mat(input1'))/1.1 SEM(cell2mat(input2'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
for i = 1:length(input1)
    cnt = cnt+1;
    plot([1 2], [nanmean(input1{i}) nanmean(input2{i})]*100, '-o', 'Color',  subj_COLOR{i}, 'MarkerEdgeColor', subj_COLOR{i})
end
ylim([0 100]); ylabel(['No response % \newline (lick ' num2str(oppwinsz) 's after reward was delivered)']);
title(['WoO = ' num2str(oppwinsz)])
xticks([1 2]); xticklabels({'Saline', 'CNO'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% uncomment if there is a third dose %%%%%%%%%%%%%%%%%%%%
% input1 = Lumped_majority_subj_performance; input2 = Lumped_treatment_subj_performance0; input3 = Lumped_treatment_subj_performance;
% figure; hold on
% bar(1:3, [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input3'), 1) nanmean(cell2mat(input2'), 1)]*100);
% errorbar(1:3, [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input3'), 1) nanmean(cell2mat(input2'), 1)]*100, ...
%     [SEM(cell2mat(input1'))/1.1 SEM(cell2mat(input3'))/1.1 SEM(cell2mat(input2'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% for i = 1:length(input1)
%     plot(1:3, [nanmean(input1{i}) nanmean(input3{i}) nanmean(input2{i})]*100, '-o', 'Color',  [1 0.844 0], 'MarkerEdgeColor', [1 0.844 0])
% end
% ylim([0 100]); ylabel('Successful withholding % \newline (lick within 1.5s after reward was delivered)');
% input1 = Lumped_majority_subj_punish; input2 = Lumped_treatment_subj_punish; input3 = Lumped_treatment_subj_punish0;
% figure; hold on
% bar(1:3, [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input3'), 1) nanmean(cell2mat(input2'), 1)]*100);
% errorbar(1:3, [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input3'), 1) nanmean(cell2mat(input2'), 1)]*100, ...
%     [SEM(cell2mat(input1'))/1.1 SEM(cell2mat(input3'))/1.1 SEM(cell2mat(input2'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% for i = 1:length(input1)
%     plot(1:3, [nanmean(input1{i}) nanmean(input3{i}) nanmean(input2{i})]*100, '-o', 'Color',  [1 0.844 0], 'MarkerEdgeColor', [1 0.844 0])
% end
% ylim([0 100]); ylabel('Failed withholding %');
% input1 = Lumped_majority_subj_slowresp; input2 = Lumped_treatment_subj_slowresp; input3 = Lumped_treatment_subj_slowresp0;
% figure; hold on;
% bar(1:3, [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input3'), 1) nanmean(cell2mat(input2'), 1)]*100);
% errorbar(1:3, [nanmean(cell2mat(input1'), 1) nanmean(cell2mat(input3'), 1) nanmean(cell2mat(input2'), 1)]*100, ...
%     [SEM(cell2mat(input1'))/1.1 SEM(cell2mat(input3'))/1.1 SEM(cell2mat(input2'))/1.1]*100, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% for i = 1:length(input1)
%     plot(1:3, [nanmean(input1{i}) nanmean(input3{i}) nanmean(input2{i})]*100, '-o', 'Color',  [1 0.844 0], 'MarkerEdgeColor', [1 0.844 0])
% end
% ylim([0 100]); ylabel('No response % \newline (lick 1.5s after reward was delivered)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
Lumped_majority_subj_RT_1 = Lumped_majority_subj_RT; 
Lumped_treatment_subj_RT_1 = Lumped_treatment_subj_RT;

for i = 1:length(Lumped_majority_subj_RT_1)
    Lumped_majority_subj_RT_1{i}(find(Lumped_majority_subj_RT_1{i}(:, 3)<0), :) = [];
end 
Lumped_majority_subj_RT_1rt = cell2mat(Lumped_majority_subj_RT_1'); Lumped_majority_subj_RT_1rt = Lumped_majority_subj_RT_1rt(:, 3);
for i = 1:length(Lumped_treatment_subj_RT_1)
    Lumped_treatment_subj_RT_1{i}(find(Lumped_treatment_subj_RT_1{i}(:, 3)<0), :) = [];
end 
Lumped_treatment_subj_RT_1rt = cell2mat(Lumped_treatment_subj_RT_1'); Lumped_treatment_subj_RT_1rt = Lumped_treatment_subj_RT_1rt(:, 3);
%%
figure; hold on
edges_RT = -1.1:0.1:10;
N = histcounts(cell2mat(Lumped_subj_RT'), edges_RT, 'Normalization', 'probability');
plot(edges_RT(2:end), N)
xlabel('Reaction time (s)')
%%
MEAN_RT_nocno = mean(Lumped_majority_subj_RT_1rt, 1); SEM_RT_nocno = SEM(cell2mat(Lumped_majority_subj_sessionRT'));
MEAN_subj_RT_nocno = nanmean(cell2mat(cellfun(@(x) nanmean(x(:, 3)), Lumped_majority_subj_RT_1, 'UniformOutput',false))); 

MEAN_RT_cno = mean(Lumped_treatment_subj_RT_1rt, 1); SEM_RT_cno = SEM(cell2mat(Lumped_treatment_subj_sessionRT'));
MEAN_subj_RT_cno = nanmean(cell2mat(cellfun(@(x) nanmean(x(:, 3)), Lumped_treatment_subj_RT_1, 'UniformOutput',false))); 
%%
figure; hold on
bar(1:2, [MEAN_subj_RT_nocno MEAN_subj_RT_cno]);
errorbar(1:2, [MEAN_RT_nocno MEAN_subj_RT_cno], [SEM_RT_nocno SEM_RT_cno]/1.3, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
RT_nocno_plt = cell2mat(Lumped_majority_subj_sessionRT');
scatter(ones(length(RT_nocno_plt), 1), RT_nocno_plt, 24, 'o', 'MarkerEdgeColor', [1 0.844 0], 'MarkerFaceColor', 'none')
RT_cno_plt = cell2mat(Lumped_treatment_subj_sessionRT');
scatter(2*ones(length(RT_cno_plt), 1), RT_cno_plt, 24, 'o', 'MarkerEdgeColor', [1 0.844 0], 'MarkerFaceColor', 'none')
for i = 1:length(Lumped_majority_subj_RT_1)
    plot([1 2], [nanmean(Lumped_majority_subj_RT_1{i}(:, 3)) nanmean(Lumped_treatment_subj_RT_1{i}(:, 3))], '-o', 'Color',  subj_COLOR{i}, 'MarkerEdgeColor', subj_COLOR{i})
end
ylabel('Reaction time (s)')
xticks([1 2]); xticklabels({'Saline', 'CNO'})