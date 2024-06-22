%%
% This script plots:
% - Tone duration-sorted performance (including baseline-normalized results)
% - Tone duration-sorted reaction time
% - Tone duration distribution (both data and simulation)
% - CNO and Saline session results can be generated 

%% Saline (or WT) running through
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 0, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 0, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 666, 0);
loopIDX = 1:30;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
OFCIDX = 1:length(Behavior_files);
len = length(OFCIDX);

num_duration_bins = 3; % number of bins used to split the [5s, 12s] tone duration
duration_bins = 5:(12-5)/num_duration_bins:12; % edges of bins
Lumped_DurationSuccess = cell(1, num_duration_bins); % binned average success rate  
Lumped_DurationRT = cell(1, num_duration_bins);         % binned average reaction time
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
OPPRT = {'RT05', 'RT075', 'RT1', 'RT15', 'RT2'};
option = 3;
cnt = 0;
for i = OFCIDX(:).'
    cnt = cnt+1;
    mx = load(MetaDataX_files{i});
    for j = 1:num_duration_bins
        MX = mx.MetaDataX(find(mx.MetaDataX.Tone_Duration>=duration_bins(j)&mx.MetaDataX.Tone_Duration<duration_bins(j+1)), :); 
        Lumped_DurationRT{j}(end+1, 1) = nanmean(MX.(OPPRT{option})(find(~isnan(MX.(OPPRT{option})))));
        Lumped_DurationSuccess{j}(end+1, 1) = length(find(MX.(OPPTWIN{option})==1))/height(MX);
    end
end

% SIM_DurationSuccess = cell(1, num_duration_bins);
% SIM_ToneDuration = cell(1, length(loopIDX)*7);
% freewin_sz = 7;
% subj_ILI = ExtractFreeILI(loopIDX, e1.MetaData, OFCIDX, ANIMAL_IDs, ANIMAL_VARs, MetaData_files, freewin_sz);
% SIM_subj_events = RI_Simulation(loopIDX, 7);
% cnt = 0;
% for i = 1:length(subj_ILI)
%     allILI = subj_ILI{i};
%     sessions_performance = NaN(length(SIM_subj_events{i}), 3);
%     for j = 1:length(SIM_subj_events{i})
%         cnt = cnt+1;
%         rng('shuffle')
%         allILI = allILI(randperm(length(allILI)));
%         alllicks = zeros(1, length(allILI)+1);
%         for n = 2:length(alllicks)
%             alllicks(n) = alllicks(n-1)+allILI(n-1);
%         end
%         endLicks = alllicks(end);
%         endEventsIDX = CrossSampling(SIM_subj_events{i}{j}, endLicks);
%         num_trials = [0 0 0];
%         num_Punish = [0 0 0];
%         for tt = 2:2:endEventsIDX
%             duration = SIM_subj_events{i}{j}(tt+1)-SIM_subj_events{i}{j}(tt);
%             SIM_ToneDuration{cnt}(end+1, 1) = duration;
%             if duration>=5 && duration<7.5
%                 if ~isempty(find(alllicks>=SIM_subj_events{i}{j}(tt)&alllicks<=SIM_subj_events{i}{j}(tt+1))==1)
%                     num_Punish(1) = num_Punish(1)+1;
%                 end
%                 num_trials(1) = num_trials(1)+1;
%             elseif duration>=7.5 && duration<10
%                 if ~isempty(find(alllicks>=SIM_subj_events{i}{j}(tt)&alllicks<=SIM_subj_events{i}{j}(tt+1))==1)
%                     num_Punish(2) = num_Punish(2)+1;
%                 end
%                 num_trials(2) = num_trials(2)+1;
%             else
%                 if ~isempty(find(alllicks>=SIM_subj_events{i}{j}(tt)&alllicks<=SIM_subj_events{i}{j}(tt+1))==1)
%                     num_Punish(3) = num_Punish(3)+1;
%                 end
%                 num_trials(3) = num_trials(3)+1;
%             end
%         end
%         sessions_performance(j, :) = 1-num_Punish./num_trials;
%     end
%     SIM_DurationSuccess{i} = sessions_performance;
% end


%% Saline running through and plot average successful withholding rate (responded before WoO) binned by tone durations
MEAN1 = NaN(length(subjIDX), num_duration_bins);
sem1 = [];
for i = 1:length(subjIDX)
    for j = 1:num_duration_bins
        MEAN1(i, j) = nanmean(Lumped_DurationSuccess{j}(subjIDX{i}))*100;
    end
end
for i = 1:num_duration_bins
    sem1(end+1) = SEM(MEAN1(:, i));
end
%%
figure;hold on
% subplot(1, 3, 1); hold on
bar(1:num_duration_bins, mean(MEAN1, 1))
for i = 1:length(subjIDX)
    plot(1:num_duration_bins, MEAN1(i, :), '-ok')
end
errorbar(1:num_duration_bins, mean(MEAN1, 1), sem1, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Success rate % WoO=1s'); ylim([0 100])
%%%%%%% plot inset panel for simulated baselines
% MEAN0 = [];
% sem0 = [];
% for i = 1:length(subjIDX)
%     MEAN0(end+1, :) = nanmean(SIM_DurationSuccess{i}, 1)*100;
% end
% for i = 1:num_duration_bins
%     sem0(end+1) = SEM(MEAN0(:, i));
% end
% axes('Position', [0.7 0.7 0.2 0.2]); hold on
% bar(1:3, mean(MEAN0))
% errorbar(1:3, mean(MEAN0), sem0, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

%%%%%%% plot tone duration-sorted successful withholding rate normalized by the simulated baseline 
% MEAN1_n = [];
% sem_n = [];
% for i = 1:length(subjIDX)
%     MEAN1_n(end+1, :) = CalcPercent(MEAN1(i, :), nanmean(SIM_DurationSuccess{i}, 1)*100);
% end
% for i = 1:num_duration_bins
%     sem_n(end+1) = SEM(MEAN1_n(:, i));
% end
% figure;hold on
% % subplot(1, 3, 1); hold on
% bar(1:num_duration_bins, mean(MEAN1_n, 1))
% for i = 1:length(subjIDX)
%     plot(1:num_duration_bins, MEAN1_n(i, :), '-ok')
% end
% errorbar(1:num_duration_bins, mean(MEAN1_n, 1), sem_n, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Performance normalized to chance level (%)'); ylim([-100 1150])
%%

%% Plot average reaction time (<WoO) binned by tone durations
figure; hold on
MEAN1 = [];
sem = [];
for i = 1:length(subjIDX)
    MEAN1(end+1, :) = [nanmean(Lumped_DurationRT{1}(subjIDX{i})) nanmean(Lumped_DurationRT{2}(subjIDX{i})) nanmean(Lumped_DurationRT{3}(subjIDX{i}))];
%     MEAN1(end+1, :) = nanmean(Lumped_DurationRT{1}(subjIDX{i}));
end
for i = 1:num_duration_bins
    sem(end+1) = SEM(MEAN1(:, i));
end
% subplot(1, 3, 1); hold on
bar(1:num_duration_bins, mean(MEAN1, 1))
for i = 1:length(subjIDX)
    plot(1:num_duration_bins, MEAN1(i, :), '-ok')
end
errorbar(1:num_duration_bins, mean(MEAN1, 1), sem, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Average Reaction Time (s)');
%% CNO running through
[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
num_duration_bins = 1; % number of bins used to split the [5s, 12s] tone duration
duration_bins = 5:(12-5)/num_duration_bins:12; % edges of bins
Lumped_DurationSuccess = cell(1, num_duration_bins); % binned average success rate  
Lumped_DurationRT = cell(1, num_duration_bins);         % binned average reaction time
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
OPPRT = {'RT05', 'RT075', 'RT1', 'RT15', 'RT2'};
option = 3;
cnt = 0;
for i = OFCIDX(:).'
    cnt = cnt+1;
    mx = load(MetaDataX_files{i});
    for j = 1:num_duration_bins
        MX = mx.MetaDataX(find(mx.MetaDataX.Tone_Duration>=duration_bins(j)&mx.MetaDataX.Tone_Duration<duration_bins(j+1)), :); 
        Lumped_DurationRT{j}(end+1, 1) = nanmean(MX.(OPPRT{option})(find(~isnan(MX.(OPPRT{option})))));
        Lumped_DurationSuccess{j}(end+1, 1) = length(find(MX.(OPPTWIN{option})==1))/height(MX);
    end
end
%% 
MEAN2 = NaN(length(subjIDX), num_duration_bins);
sem2 = [];
for i = 1:length(subjIDX)
    for j = 1:num_duration_bins
        MEAN2(i, j) = nanmean(Lumped_DurationSuccess{j}(subjIDX{i}))*100;
    end
end
for i = 1:num_duration_bins
    sem2(end+1) = SEM(MEAN2(:, i));
end
%%
subplot(1, 3, 2); hold on
bar(1:num_duration_bins, mean(MEAN2, 1))
for i = 1:length(subjIDX)
    plot(1:num_duration_bins, MEAN2(i, :), '-ok')
end
errorbar(1:num_duration_bins, mean(MEAN2, 1), sem, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Average % Successful Trials (reaction time < 1.5s)'); ylim([0 100])
%% Plot average reaction time (<WoO) binned by tone durations
MEAN2 = [];
sem = [];
for i = 1:length(subjIDX)
%     MEAN2(end+1, :) = [nanmean(Lumped_DurationRT{1}(subjIDX{i})) nanmean(Lumped_DurationRT{2}(subjIDX{i})) nanmean(Lumped_DurationRT{3}(subjIDX{i}))];
    MEAN2(end+1, :) = nanmean(Lumped_DurationRT{1}(subjIDX{i}));
end
for i = 1:num_duration_bins
    sem(end+1) = SEM(MEAN2(:, i));
end
subplot(1, 3, 2); hold on
bar(1:num_duration_bins, mean(MEAN2, 1))
for i = 1:length(subjIDX)
    plot(1:num_duration_bins, MEAN2(i, :), '-ok')
end
errorbar(1:num_duration_bins, mean(MEAN2, 1), sem, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Average Reaction Time (s)');
%% Plot successful withholding rate between Saline and CNO group binned by tone durations
figure; hold on
for i = 1:num_duration_bins
    bar(3*i-2, mean(MEAN1(:, i)))
    bar(3*i-1, mean(MEAN2(:, i)))
    for j = 1:size(MEAN1, 1)
        plot(3*i-2:3*i-1, [MEAN1(j, i) MEAN2(j, i)], '-ok')
    end
    errorbar(3*i-2:3*i-1, [mean(MEAN1(:, i)) mean(MEAN2(:, i))], [sem1(i) sem2(i)], 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
end
ylabel('Success rate (%)')
xlabel('Inhibition Tone Duration (s)')






%% Plot ratio of average successful withholding rate between Saline and CNO group binned by tone durations
MEAN_ratio = MEAN2./MEAN1;
sem = [];
for i = 1:num_duration_bins
    sem(end+1) = SEM(MEAN_ratio(:, i));
end
subplot(1, 3, 3); hold on
bar(1:num_duration_bins, mean(MEAN_ratio, 1))
for i = 1:length(subjIDX)
    plot(1:num_duration_bins, MEAN_ratio(i, :), '-ok')
end
errorbar(1:num_duration_bins, mean(MEAN_ratio, 1), sem, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
%% Plot ratio of average reaction time Saline vs CNO binned by tone durations
MEAN_RT1 = mean(MEAN1, 2);
MEAN_RT2 = mean(MEAN2, 2);
sem1 = SEM(MEAN_RT1);
sem2 = SEM(MEAN_RT2);
% subplot(1, 3, 3); hold on
figure; hold on
bar(1, mean(MEAN_RT1))
bar(2, mean(MEAN_RT2))
for i = 1:length(subjIDX)
    plot(1:2, [MEAN_RT1(i) MEAN_RT2(i)], '-ok')
end
errorbar(1:2, [mean(MEAN_RT1) mean(MEAN_RT2)], [sem1 sem2], 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

%% Tone duration distribution of simulation
num_bins = 20;
MEAN_ToneDuration_distr = [];
SEM_ToneDuration_distr = [];
for i = 1:length(SIM_ToneDuration)
    [N,edges] = histcounts(SIM_ToneDuration{i}, num_bins, 'Normalization', 'probability');
    MEAN_ToneDuration_distr(end+1, :) = N;
end
for i = 1:num_bins
    SEM_ToneDuration_distr(end+1) = SEM(MEAN_ToneDuration_distr(:, i));
end
figure; hold on
bar(5:(12-5)/num_bins:12-(12-5)/num_bins, mean(MEAN_ToneDuration_distr, 1))
errorbar(5:(12-5)/num_bins:12-(12-5)/num_bins, mean(MEAN_ToneDuration_distr, 1), SEM_ToneDuration_distr, 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Fraction'); xlabel('Tone Duration (s)')
