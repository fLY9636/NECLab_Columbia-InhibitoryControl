%% Create column vectors containing all sessions' NE and Ach traces 
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 1);
%****************** for temporary test **************************
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 7777, 0);
loopIDX = 1:2;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
OFCIDX = 1:length(Behavior_files);
len = length(OFCIDX);
e1.MetaData = MetaData_GRABmute;
%****************************************************************

IDX = 1500*120;
Lumped_NE = NaN(IDX, len);
Lumped_Ach = NaN(IDX, len);
for i = OFCIDX(:).'
    Lumped_NE(:, i) = zscore(Filter(e1.MetaData.NE_470{i}(1:IDX), 120, 4, [0.1 3.5], 'bandpass'));
%     Lumped_NE(:, i) = Filter(e1.MetaData.NE_corr{i}(1:IDX), 120, 4, [0.1 3.5], 'bandpass');
    Lumped_Ach(:, i) = zscore(Filter(e1.MetaData.Ach_470{i}(1:IDX), 120, 4, [0.1 3.5], 'bandpass'));
%     Lumped_Ach(:, i) = Filter(e1.MetaData.Ach_corr{i}(1:IDX), 120, 4, [0.1 3.5], 'bandpass');
end

%% Parameters initialization
% See detailed documentation on http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/continuous/cohgramc.html
% movingwin = [median(EngageDuration) median(EngageDuration)/2];
movingwin = [11 6]; % 
params.pad = 0;
params.tapers = [5 8];
params.Fs = 120;
params.fpass = [0 60];
% params.err = [0 0.05];
params.trialave = 0;
%% Coherence calculation
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 8);
figure; hold on
[C_Lumped, phi, S12, S1, S2, t, f1] = cohgramc(Lumped_NE(:, OFCIDX), Lumped_Ach(:, OFCIDX), movingwin, params); % get the lumped sessions coherence traces
% uncomment plot the averaged result only when params.trialave=1
% [C, phi, S12, S1, S2, t, f] = cohgramc(Lumped_NE(:, OFCIDX), Lumped_Ach(:, OFCIDX), movingwin, params);
% plot(f, nanmean(C, 1))
% xlabel('Frequency'); ylabel('Coherence');
% set(gca, 'XScale', 'Log')


figure; hold on
subjC_set = [];
cnt = 0;
% plot averaged coherence across subjects and individual coherence 
for subj = 1:length(subjIDX)
    sessionmeanC_set = [];
    cnt = cnt+1;
    for i = subjIDX{subj}(:).'
        cnt = cnt+1;
        sessionmeanC_set(end+1, :) = nanmean(reshape(C_Lumped(:, :, i), size(C_Lumped, 1), size(C_Lumped, 2)), 1);
    end
    subjC_set(end+1, :) = nanmean(sessionmeanC_set, 1);
%     plot(f1, nanmean(sessionmeanC_set, 1), 'k', 'LineWidth', 0.5)
end
idx = find(f1>=0.05);
ShadedPlot(f1(idx), nanmean(subjC_set(:, idx), 1), [0 0 1], 2, SEM(subjC_set(:, idx)), [0.8 0.9 0.9])
plot(f1, nanmean(subjC_set, 1), 'b', 'LineWidth', 1)
xlim([0.05 3.5])

pltIDX = find(f1>=0.05);
h = patch([f1(pltIDX) fliplr(f1(pltIDX))], [QL fliplr(QU)], [0, .6, .77], 'FaceAlpha',0.5, 'EdgeColor','none')