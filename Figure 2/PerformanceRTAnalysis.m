% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
%%%%%%%%%%% for temporary folder run-through ************
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
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task performance and average reaction time of each session
Lumped_Performance = NaN(length(OFCIDX), 1);
Lumped_RT = NaN(length(OFCIDX), 1);
LumpedTrial_RT = [];
for i = OFCIDX(:).'
    mx = load(MetaDataX_files{i});
    Lumped_Performance(i) = length(find(mx.MetaDataX.Outcome1==1))/height(mx.MetaDataX);
    Lumped_RT(i) = nanmean(mx.MetaDataX.RT1);
    LumpedTrial_RT = [LumpedTrial_RT; mx.MetaDataX.RT1(~isnan(mx.MetaDataX.RT1))];
end
%%
num_bins = 4;
subj_Performance_RT_pair = NaN(length(subjIDX), num_bins);
for i = 1:length(subjIDX)
    subj_RT = Lumped_RT(subjIDX{i});
    subj_Performance = Lumped_Performance(subjIDX{i});
    [MeanINBin, SEMINBin, bin_ctrs, bin_IDXs] = BinnedQTMean1d(subj_RT, subj_Performance, num_bins);
    subj_Performance_RT_pair(i, :) = MeanINBin;
end
figure; hold on
plot(bin_ctrs, nanmean(subj_Performance_RT_pair, 1))
for i = 1:length(subjIDX)
    plot(bin_ctrs, subj_Performance_RT_pair(i, :), '-ok')
end
errorbar(bin_ctrs, nanmean(subj_Performance_RT_pair, 1), SEM(subj_Performance_RT_pair), 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
xlabel('Percentile of individual reaction time %')
ylabel('Average success rate')
%%
num_bins = 3;
norm_RT = UnitNormalization(Lumped_RT);
subj_RT = Lumped_RT;
subj_Performance = Lumped_Performance;
[MeanINBin, SEMINBin, bin_ctrs, bin_IDXs] = BinnedMean1d(subj_RT, subj_Performance, num_bins);
figure; hold on
plot(bin_ctrs, MeanINBin)
errorbar(bin_ctrs, MeanINBin, SEMINBin, 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

%%
figure; hold on
scatter(Lumped_RT, Lumped_Performance)
%% Plot reaction time distribution - supplementary figure 2e
edges = 0:0.02:1;
% figure; hold on
% for i = 1:length(subjIDX)
%     subplot(4, 5, i); hold on
%     histogram(LumpedTrial_RT(subjIDX{i}), edges)
%     title(ANIMAL_IDs{i})
% end
figure; hold on
% subplot(2, 3, 3)
% histogram(LumpedTrial_RT, edges)
% subplot(2, 3, 6)
histogram(LumpedTrial_RT, edges, 'Normalization', 'probability')
%%
% gmm = fitgmdist(LumpedTrial_RT, 2, );
% idx = cluster(gmm, LumpedTrial_RT);
% cluster1 = LumpedTrial_RT(idx == 1);
% cluster2 = LumpedTrial_RT(idx == 2);
%%
edges = 0:0.01:1;
figure; hold on
histogram(cluster1, edges)
histogram(cluster2, edges)
