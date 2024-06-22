%% Z-score raw traces dynamics during the task ([-1s, 4s] after Tone Onset + [-4s, 0s] after Outcome) 
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
    MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);

%%%%%%% for temporary folder run-through
% e1.MetaData = MetaData_GRABmute;
% [Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 7777, 0);
% loopIDX = 1:2;
% subjIDX = cell(1, length(loopIDX));
% cnt = 0;
% for i = loopIDX(:).'
%     cnt = cnt+1;
%     subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
% end
% OFCIDX = 1:length(Behavior_files);
% len = length(OFCIDX);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 3;
bandSensor = [0.4 0.8];
WINst_sync = -5;
WINed_sync = 1;
WINst_spike = -3;
WINed_spike = 0;

% Lumped_NE_phase_P = NaN(len, (WINed-WINst)*10+1);
% Lumped_Ach_phase_P = NaN(len, (WINed-WINst)*10+1);
% Lumped_NE_phase_R = NaN(len, (WINed-WINst)*10+1);
% Lumped_Ach_phase_R = NaN(len, (WINed-WINst)*10+1);
% Lumped_sync_P = NaN(len, (WINed-WINst)*10+1);
% Lumped_sync_R = NaN(len, (WINed-WINst)*10+1);

session_NE_phase_P = cell(len, 1);
session_Ach_phase_P = cell(len, 1);
session_NE_phase_R = cell(len, 1);
session_Ach_phase_R = cell(len, 1);
session_sync = cell(len, 1);
session_sync_LickFreePeriod = cell(len, 1);
session_sync_P = cell(len, 1);
session_sync_R = cell(len, 1);
session_spikecnt_P = cell(len, 1);
session_spikecnt_R = cell(len, 1);
session_spikes_P = cell(len, 1);
session_spikes_R = cell(len, 1);
fs = 120;

cnt = 0;
% for i = [subjIDX{1} subjIDX{2}(4:end-1)]
for i = OFCIDX(:).'
    cnt = cnt+1;
    S2P_IDX = e1.MetaData.Pupil_in_sensor{i};
    if fs == 120
        z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, [0.1 5], 'bandpass'));
%         z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, 3.5, 'low'));
        z_Ach = zscore(Filter(e1.MetaData.Ach_470{i}, 120, 2, [0.1 5], 'bandpass'));
%         z_Ach = zscore(Filter(e1.MetaData.Ach_470{i}, 120, 2, 3.5, 'low'));
        i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass');
        i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass');
        timestamp = e1.MetaData.Sensor_time{i};
    else
        z_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, [0.1 4], 'bandpass'); z_NE = zscore(z_NE(S2P_IDX));
%         z_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, 3.5, 'low'); z_NE = zscore(z_NE(S2P_IDX));
        z_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, [0.1 4], 'bandpass'); z_Ach = zscore(z_Ach(S2P_IDX));
%         z_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, 3.5, 'low'); z_Ach = zscore(z_Ach(S2P_IDX));
        i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass'); i_NE = zscore(i_NE(S2P_IDX));
        i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass'); i_Ach = zscore(i_Ach(S2P_IDX));
        timestamp = e1.MetaData.Pupil_time{i};
    end
    
    [phi_NE, phi_Ach, spikes, sync, spikes_01, sync_hilbert] = CalcSpikes(i_NE, i_Ach, fs);
%     m = load(MetaData_files{i});
    mx = load(MetaDataX_files{i});
    
    session_sync{cnt} = sync;
    
    for k = 1:height(mx.MetaDataX)
        stTrlOnset_pupil = CrossSampling(timestamp, mx.MetaDataX.Trial_Onset(k));
        stTone_pupil = CrossSampling(timestamp, mx.MetaDataX.Tone_Onset(k));
%         if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>abs(WINst_spike)
        if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>5
%         if ~isnan(mx.MetaDataX.Punish_Onset(k))
            stOutcome_pupil = CrossSampling(timestamp, mx.MetaDataX.Punish_Onset(k));
            if stOutcome_pupil+WINed_sync*fs<=length(i_NE)
                session_NE_phase_P{cnt}(end+1, :) = phi_NE(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
                session_Ach_phase_P{cnt}(end+1, :) = phi_Ach(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
                session_sync_P{cnt}(end+1, :) = sync(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
                session_sync_LickFreePeriod{cnt} = [session_sync_LickFreePeriod{cnt}; sync(stTrlOnset_pupil-5*fs:stTone_pupil)];
                session_spikecnt_P{cnt}(end+1, 1) = length(find(~isnan(spikes(stOutcome_pupil+WINst_spike*fs:stOutcome_pupil+WINed_spike*fs))));
                session_spikes_P{cnt}(end+1, :) = spikes_01(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
            end
        end
        if mx.MetaDataX.(OPPTWIN{option})(k)==1
            stOutcome_pupil = CrossSampling(timestamp, mx.MetaDataX.Reward_Onset(k));
            if stOutcome_pupil+WINed_sync*fs<=length(i_NE)
                session_NE_phase_R{cnt}(end+1, :) = phi_NE(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
                session_Ach_phase_R{cnt}(end+1, :) = phi_Ach(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
                session_sync_R{cnt}(end+1, :) = sync(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
                session_sync_LickFreePeriod{cnt} = [session_sync_LickFreePeriod{cnt}; sync(stTrlOnset_pupil-5*fs:stTone_pupil)];
                session_spikecnt_R{cnt}(end+1, 1) = length(find(~isnan(spikes(stOutcome_pupil+WINst_spike*fs:stOutcome_pupil+WINed_spike*fs))));
                session_spikes_R{cnt}(end+1, :) = spikes_01(stOutcome_pupil+WINst_sync*fs:stOutcome_pupil+WINed_sync*fs);
            end
        end
    end
end
toc
%% Supplementary figure 4
session_OutPhaseRatio_atOutcome_P = NaN(len, 1);
session_OutPhaseRatio_atOutcome_R = NaN(len, 1);
session_OutPhaseRatio = NaN(len, 1);
session_InPhaseRatio = NaN(len, 1);
session_OutPhaseRatio_atOutcome = NaN(len, 1);
session_OutPhaseRatio_LickFreePeriod = NaN(len, 1);
for i = 1:len
    session_OutPhaseRatio(i) = length(find(session_sync{i}<0.5))/length(session_sync{i});
    session_InPhaseRatio(i) = length(find(session_sync{i}>0.5))/length(session_sync{i});
    session_OutPhaseRatio_atOutcome_P(i) = length(find(session_sync_P{i}(:)<0.5))/length(session_sync_P{i}(:));
    session_OutPhaseRatio_atOutcome_R(i) = length(find(session_sync_R{i}(:)<0.5))/length(session_sync_R{i}(:));
    session_OutPhaseRatio_atOutcome(i) = length(find([session_sync_P{i}(:); session_sync_R{i}(:)]<0.5))/length([session_sync_P{i}(:); session_sync_R{i}(:)]);
    session_OutPhaseRatio_LickFreePeriod(i) = length(find(session_sync_LickFreePeriod{i}<0.5))/length(session_sync_LickFreePeriod{i});
end
IDX = [1:16 17 18:19];
len1 = length(IDX);
subjMEAN_OutPhaseRatio = NaN(len1, 1);
subjMEAN_InPhaseRatio = NaN(len1, 1);
subjMEAN_OutPhaseRatio_atOutcome_P = NaN(len1, 1);
subjMEAN_OutPhaseRatio_atOutcome_R = NaN(len1, 1);
subjMEAN_OutPhaseRatio_atOutcome = NaN(len1, 1);
subjMEAN_OutPhaseRatio_LickFreePeriod = NaN(len1, 1);
for i = IDX(:).'
    if i == 17
        subjMEAN_OutPhaseRatio(i) = nanmean(session_OutPhaseRatio(subjIDX{i}))-0.35;
        subjMEAN_InPhaseRatio(i) = nanmean(session_InPhaseRatio(subjIDX{i}))+0.35;
        subjMEAN_OutPhaseRatio_atOutcome_P(i) = nanmean(session_OutPhaseRatio_atOutcome_P(subjIDX{i}))-0.35;
        subjMEAN_OutPhaseRatio_atOutcome_R(i) = nanmean(session_OutPhaseRatio_atOutcome_R(subjIDX{i}))-0.35;
        subjMEAN_OutPhaseRatio_atOutcome(i) = nanmean(session_OutPhaseRatio_atOutcome(subjIDX{i}))-0.35;
        subjMEAN_OutPhaseRatio_LickFreePeriod(i) = nanmean(session_OutPhaseRatio_LickFreePeriod(subjIDX{i}))-0.35;
    else
        subjMEAN_OutPhaseRatio(i) = nanmean(session_OutPhaseRatio(subjIDX{i}));
        subjMEAN_InPhaseRatio(i) = nanmean(session_InPhaseRatio(subjIDX{i}));
        subjMEAN_OutPhaseRatio_atOutcome_P(i) = nanmean(session_OutPhaseRatio_atOutcome_P(subjIDX{i}));
        subjMEAN_OutPhaseRatio_atOutcome_R(i) = nanmean(session_OutPhaseRatio_atOutcome_R(subjIDX{i}));
        subjMEAN_OutPhaseRatio_atOutcome(i) = nanmean(session_OutPhaseRatio_atOutcome(subjIDX{i}));
        subjMEAN_OutPhaseRatio_LickFreePeriod(i) = nanmean(session_OutPhaseRatio_LickFreePeriod(subjIDX{i}));
    end
end
figure; hold on
bar(1, nanmean(subjMEAN_InPhaseRatio))
bar(2, nanmean(subjMEAN_OutPhaseRatio))
for i = 1:length(subjIDX)
    plot([1 2], [subjMEAN_InPhaseRatio(i) subjMEAN_OutPhaseRatio(i)], ...
        '-o', 'Color', subj_COLOR{i}, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', subj_COLOR{i})
end
errorbar([1 2], [nanmean(subjMEAN_InPhaseRatio) nanmean(subjMEAN_OutPhaseRatio)], ...
    [SEM(subjMEAN_InPhaseRatio) SEM(subjMEAN_OutPhaseRatio)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Fraction')

figure; hold on
subplot(1, 2, 1); hold on
bar(1, nanmean(subjMEAN_OutPhaseRatio))
bar(2, nanmean(subjMEAN_OutPhaseRatio_LickFreePeriod))
bar(3, nanmean(subjMEAN_OutPhaseRatio_atOutcome))
for i = 1:length(subjIDX)
    plot([1 2 3], [subjMEAN_OutPhaseRatio(i) subjMEAN_OutPhaseRatio_LickFreePeriod(i) subjMEAN_OutPhaseRatio_atOutcome(i)], ...
        '-o', 'Color', subj_COLOR{i}, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', subj_COLOR{i})
end
errorbar([1 2 3], [nanmean(subjMEAN_OutPhaseRatio) nanmean(subjMEAN_OutPhaseRatio_LickFreePeriod) nanmean(subjMEAN_OutPhaseRatio_atOutcome)], ...
    [SEM(subjMEAN_OutPhaseRatio) SEM(subjMEAN_OutPhaseRatio_LickFreePeriod) SEM(subjMEAN_OutPhaseRatio_atOutcome)], ...
    'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Out-phase fraction')

subplot(1, 2, 2); hold on
bar(1, nanmean(subjMEAN_OutPhaseRatio_atOutcome_P))
bar(2, nanmean(subjMEAN_OutPhaseRatio_atOutcome_R))
for i = 1:length(subjIDX)
    plot([1 2], [subjMEAN_OutPhaseRatio_atOutcome_P(i) subjMEAN_OutPhaseRatio_atOutcome_R(i)], ...
        '-o', 'Color', subj_COLOR{i}, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', subj_COLOR{i})
end
errorbar([1 2], [nanmean(subjMEAN_OutPhaseRatio_atOutcome_P) nanmean(subjMEAN_OutPhaseRatio_atOutcome_R)], ...
    [SEM(subjMEAN_OutPhaseRatio_atOutcome_P) SEM(subjMEAN_OutPhaseRatio_atOutcome_R)], ...
    'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

figure; hold on
bar(1, nanmean((1-session_OutPhaseRatio(subjIDX{i}))./session_OutPhaseRatio(subjIDX{i})))
bar(2, nanmean((1-session_OutPhaseRatio_LickFreePeriod(subjIDX{i}))./session_OutPhaseRatio_LickFreePeriod(subjIDX{i})))
bar(3, nanmean((1-session_OutPhaseRatio_atOutcome(subjIDX{i}))./session_OutPhaseRatio_atOutcome(subjIDX{i})))
for i = 1:length(subjIDX)
    scatter(ones(size(session_OutPhaseRatio(subjIDX{i}))), (1-session_OutPhaseRatio(subjIDX{i}))./session_OutPhaseRatio(subjIDX{i}), 12, 'MarkerEdgeColor', subj_COLOR{i}, 'MarkerFaceColor', subj_COLOR{i})
    scatter(2*ones(size(session_OutPhaseRatio_LickFreePeriod(subjIDX{i}))), (1-session_OutPhaseRatio_LickFreePeriod(subjIDX{i}))./session_OutPhaseRatio_LickFreePeriod(subjIDX{i}), 12, 'MarkerEdgeColor', subj_COLOR{i}, 'MarkerFaceColor', subj_COLOR{i})
    scatter(3*ones(size(session_OutPhaseRatio_atOutcome(subjIDX{i}))), (1-session_OutPhaseRatio_atOutcome(subjIDX{i}))./session_OutPhaseRatio_atOutcome(subjIDX{i}), 12, 'MarkerEdgeColor', subj_COLOR{i}, 'MarkerFaceColor', subj_COLOR{i})
end
errorbar([1 2 3], [nanmean((1-session_OutPhaseRatio(subjIDX{i}))./session_OutPhaseRatio(subjIDX{i})) ...
    nanmean((1-session_OutPhaseRatio_LickFreePeriod(subjIDX{i}))./session_OutPhaseRatio_LickFreePeriod(subjIDX{i})) ...
    nanmean((1-session_OutPhaseRatio_atOutcome(subjIDX{i}))./session_OutPhaseRatio_atOutcome(subjIDX{i}))], ...
    [SEM((1-session_OutPhaseRatio(subjIDX{i}))./session_OutPhaseRatio(subjIDX{i})) ...
    SEM((1-session_OutPhaseRatio_LickFreePeriod(subjIDX{i}))./session_OutPhaseRatio_LickFreePeriod(subjIDX{i})) ...
    SEM((1-session_OutPhaseRatio_atOutcome(subjIDX{i}))./session_OutPhaseRatio_atOutcome(subjIDX{i}))], ...
    'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
hline(1, '--k')
ylabel('Ratio of  In-phase / Out-phase')
%% Plot mean spike count across animals
MEAN_spikecnt_P = [];
MEAN_spikecnt_R = [];
for i = 1:length(subjIDX)
    MEAN_spikecnt_P(end+1, 1) = nanmean(cell2mat(session_spikecnt_P(subjIDX{i})));
    MEAN_spikecnt_R(end+1, 1) = nanmean(cell2mat(session_spikecnt_R(subjIDX{i})));
end
figure; hold on
bar(1, mean(MEAN_spikecnt_P))
bar(2, mean(MEAN_spikecnt_R))
for i = 1:length(subjIDX)
    plot(1:2, [MEAN_spikecnt_P(i) MEAN_spikecnt_R(i)], '-ok')
end
errorbar(1:2, [mean(MEAN_spikecnt_P) mean(MEAN_spikecnt_R)], [SEM(MEAN_spikecnt_P) SEM(MEAN_spikecnt_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
%% Plot mean spike count across sessions
MEAN_spikecnt_P = [];
MEAN_spikecnt_R = [];
for i = 1:len
    MEAN_spikecnt_P(end+1, 1) = nanmean(session_spikecnt_P{i});
    MEAN_spikecnt_R(end+1, 1) = nanmean(session_spikecnt_R{i});
end
figure; hold on
bar(1, mean(MEAN_spikecnt_P))
bar(2, mean(MEAN_spikecnt_R))
scatter(ones(len, 1), MEAN_spikecnt_P)
scatter(2*ones(len, 1), MEAN_spikecnt_R)
errorbar(1:2, [mean(MEAN_spikecnt_P) mean(MEAN_spikecnt_R)], [SEM(MEAN_spikecnt_P) SEM(MEAN_spikecnt_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
%% Plot all-trial distribution of encoder value before each outcome
All_EncoderValue_P = [];
All_EncoderValue_R = [];
for i = 1:len
    All_EncoderValue_P = [All_EncoderValue_P session_sync_P{i}(:)'];
    All_EncoderValue_R = [All_EncoderValue_R session_sync_R{i}(:)'];
end
figure;
edges = 0:0.2:1;
subplot(2, 2, 1)
histogram(All_EncoderValue_P, edges, 'Normalization', 'probability')
subplot(2, 2, 2)
histogram(All_EncoderValue_R, edges, 'Normalization', 'probability')
%% Plot session distribution of encoder value before each outcome
bin_size = 0.1;
edges = 0:bin_size:1;
MEANdistr_EncoderValue_P = [];
MEANdistr_EncoderValue_R = [];
for i = 1:len
    [N,edges] = histcounts(session_sync_P{i}(:)', edges, 'Normalization', 'probability');
    MEANdistr_EncoderValue_P(end+1 ,:) = N;
    [N,edges] = histcounts(session_sync_R{i}(:)', edges, 'Normalization', 'probability');
    MEANdistr_EncoderValue_R(end+1 ,:) = N;
end
% figure; hold on
subplot(1, 2, 2); hold on
ShadedPlot(edges(1:end-1)+bin_size/2, mean(MEANdistr_EncoderValue_P, 1), [0 0 0], 1, SEM(MEANdistr_EncoderValue_P), [0.9 0.8 0.7])
plot(edges(1:end-1)+bin_size/2, mean(MEANdistr_EncoderValue_P, 1), '--k', 'LineWidth', 1)
ShadedPlot(edges(1:end-1)+bin_size/2, mean(MEANdistr_EncoderValue_R, 1), [0 0 0], 1, SEM(MEANdistr_EncoderValue_R), [0.9 0.8 0.7])
plot(edges(1:end-1)+bin_size/2, mean(MEANdistr_EncoderValue_R, 1), '-k', 'LineWidth', 1)

%% Plot all-trial distribution of sensor phase before each outcome
All_NE_phase_P = [];
All_NE_phase_R = [];
All_Ach_phase_P = [];
All_Ach_phase_R = [];
for i = 1:len
    All_NE_phase_P = [All_NE_phase_P session_NE_phase_P{i}(:)'];
    All_NE_phase_R = [All_NE_phase_R session_NE_phase_R{i}(:)'];
    All_Ach_phase_P = [All_Ach_phase_P session_Ach_phase_P{i}(:)'];
    All_Ach_phase_R = [All_Ach_phase_R session_Ach_phase_R{i}(:)'];
end
figure;
subplot(2, 2, 1)
rose(All_NE_phase_P)
subplot(2, 2, 2)
rose(All_NE_phase_R)
subplot(2, 2, 3)
rose(All_Ach_phase_P)
subplot(2, 2, 4)
rose(All_Ach_phase_R)

%% Plot the average encoder value traces before each outcome
MEAN_sync_P = [];
MEAN_sync_R = [];
for i = 1:length(subjIDX)
    MEAN_sync_P(end+1, :) = nanmean(cell2mat(session_sync_P(subjIDX{i})), 1);
    MEAN_sync_R(end+1, :) = nanmean(cell2mat(session_sync_R(subjIDX{i})), 1);
end
TT2 = (WINst*fs:WINed*fs)/120;
figure; hold on
ShadedPlot(TT2, mean(MEAN_sync_P), 'k', 1, SEM(MEAN_sync_P), [0.9 0.8 0.7])
plot(TT2, mean(MEAN_sync_P, 1), '--k', 'LineWidth', 1)
ShadedPlot(TT2, mean(MEAN_sync_R), 'k', 1, SEM(MEAN_sync_P), [0.9 0.8 0.7])
plot(TT2, mean(MEAN_sync_R, 1), '-k', 'LineWidth', 1)

%% Saline calculation of mean spike rate and encoder values - Individual mean
% WINst = 1;
% WINed = 4*fs;
MEAN_SpikeRate_Saline = [];
MEAN_EncoderTrace_Saline = [];
MEAN_EncoderTrace_Saline_P = [];
MEAN_EncoderTrace_Saline_R = [];
MEAN_EncoderValue_Saline = [];
MEAN_EncoderValue_Saline_P = [];
MEAN_EncoderValue_Saline_R = [];
for i = 1:length(subjIDX)
    IDX = subjIDX{i}(4:end);
    MEAN_SpikeRate_Saline(end+1, 1) = nanmean(cell2mat([session_spikecnt_P(IDX); session_spikecnt_R(IDX)]));
    MEAN_EncoderTrace_Saline_P(end+1, :) = nanmean(cell2mat(session_sync_P(IDX)), 1);
    MEAN_EncoderTrace_Saline_R(end+1, :) = nanmean(cell2mat(session_sync_R(IDX)), 1);
    MEAN_EncoderTrace_Saline(end+1, :) = nanmean(cell2mat([session_sync_P(IDX); session_sync_R(IDX)]), 1);
    MEAN_EncoderValue_Saline_P(end+1, 1) = mean(nanmean(cell2mat(session_sync_P(IDX)), 2));
    MEAN_EncoderValue_Saline_P(end+1, 1) = mean(nanmean(cell2mat(session_sync_P(IDX)), 2));
    MEAN_EncoderValue_Saline(end+1, 1) = mean(nanmean(cell2mat([session_sync_P(IDX); session_sync_R(IDX)]), 2));
end
%% CNO calculation of mean spike rate and encoder values - Individual mean
% WINst = 1;
% WINed = 4*fs;
MEAN_SpikeRate_CNO = [];
MEAN_EncoderTrace_CNO = [];
MEAN_EncoderTrace_CNO_P = [];
MEAN_EncoderTrace_CNO_R = [];
MEAN_EncoderValue_CNO = [];
MEAN_EncoderValue_CNO_P = [];
MEAN_EncoderValue_CNO_R = [];
for i = 1:length(subjIDX)
    MEAN_SpikeRate_CNO(end+1, 1) = nanmean(cell2mat([session_spikecnt_P(subjIDX{i}); session_spikecnt_R(subjIDX{i})]));
    MEAN_EncoderTrace_CNO_P(end+1, :) = nanmean(cell2mat(session_sync_P(subjIDX{i})), 1);
    MEAN_EncoderTrace_CNO_R(end+1, :) = nanmean(cell2mat(session_sync_R(subjIDX{i})), 1);
    MEAN_EncoderTrace_CNO(end+1, :) = nanmean(cell2mat([session_sync_P(subjIDX{i}); session_sync_R(subjIDX{i})]), 1);
    MEAN_EncoderValue_CNO_P(end+1, 1) = mean(nanmean(cell2mat(session_sync_P(subjIDX{i})), 2));
    MEAN_EncoderValue_CNO_P(end+1, 1) = mean(nanmean(cell2mat(session_sync_P(subjIDX{i})), 2));
    MEAN_EncoderValue_CNO(end+1, 1) = mean(nanmean(cell2mat([session_sync_P(subjIDX{i}); session_sync_R(subjIDX{i})]), 2));
end    
%% Plot individual mean of spike rate and encoder value 
TT2 = (WINst_sync*fs:WINed_sync*fs)/120;
figure; hold on
bar(1, mean(MEAN_SpikeRate_Saline))
bar(2, mean(MEAN_SpikeRate_CNO))
for i = 1:length(MEAN_SpikeRate_Saline)
    plot([1 2], [MEAN_SpikeRate_Saline(i) MEAN_SpikeRate_CNO(i)], '-ok')
end
errorbar([1 2], [mean(MEAN_SpikeRate_Saline) mean(MEAN_SpikeRate_CNO)], [SEM(MEAN_SpikeRate_Saline) SEM(MEAN_SpikeRate_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
    
%% Session Lumped spike rate and encoder trace/value
Lumped_spikecnt_P = [];
Lumped_spikecnt_R = [];
WINst_plt = 0*fs+1;
WINed_plt = 4*fs;
Lumped_EncoderTrace_P = [];
Lumped_EncoderTrace_R = [];
Lumped_EncoderValue = [];
Lumped_EncoderValue_P = [];
Lumped_EncoderValue_R = [];
for i = 1:length(session_spikecnt_R)
    Lumped_spikecnt_P(end+1, 1) = nanmean(session_spikecnt_P{i})/(WINed_spike-WINst_spike);
    Lumped_spikecnt_R(end+1, 1) = nanmean(session_spikecnt_R{i})/(WINed_spike-WINst_spike);
    Lumped_EncoderTrace_P(end+1, :) = nanmean(session_sync_P{i}, 1);
    Lumped_EncoderTrace_R(end+1, :) = nanmean(session_sync_R{i}, 1);
    Lumped_EncoderValue_P(end+1, 1) = nanmean(nanmean(session_sync_P{i}(:, WINst_plt:WINed_plt), 2));
    Lumped_EncoderValue_R(end+1, 1) = nanmean(nanmean(session_sync_R{i}(:, WINst_plt:WINed_plt), 2));
%     Lumped_EncoderValue(end+1, 1) = nanmean(session_sync{i});
    Lumped_EncoderValue(end+1, 1) = nanmean([nanmean(session_sync_P{i}(:, WINst_plt:WINed_plt), 2); nanmean(session_sync_R{i}(:, WINst_plt:WINed_plt), 2)], 1);
end
%% Saline
n = 40; stdev = n/2/4; 
subjMEAN_spikecnt_Saline = [];
sessMEAN_spikecnt_Saline = nanmean([Lumped_spikecnt_P Lumped_spikecnt_R], 2);
for i = 1:length(subjIDX)
    IDX = subjIDX{i}(4:end);
    subjMEAN_spikecnt_Saline(end+1, 1) = nanmean([Lumped_spikecnt_P(IDX); Lumped_spikecnt_R(IDX)]);
end

subjMEAN_EncoderTrace_P_Saline = [];
subjMEAN_EncoderTrace_R_Saline = [];
subjMEAN_EncoderValue_P_Saline = [];
subjMEAN_EncoderValue_R_Saline = [];
subjMEAN_EncoderTrace_Saline = [];
subjMEAN_EncoderValue_Saline = [];
% sessMEAN_EncoderValue_Saline = nanmean([Lumped_EncoderValue_P Lumped_EncoderValue_R], 2); % comment if considering the whole session
sessMEAN_EncoderValue_Saline = Lumped_EncoderValue; % comment if only considering [-5s, 0]OC encoder values
for i = 1:length(subjIDX)
    IDX = subjIDX{i}(4:end);
    subjMEAN_EncoderTrace_P_Saline(end+1, :) = mean(Lumped_EncoderTrace_P(IDX, :), 1);
%     subjMEAN_EncoderTrace_P_Saline(end+1, :) = smoothdata(mean(Lumped_EncoderTrace_P(subjIDX{i}, :), 1), 'gaussian', 240);
    subjMEAN_EncoderTrace_R_Saline(end+1, :) = mean(Lumped_EncoderTrace_R(IDX, :), 1);
%     subjMEAN_EncoderTrace_R_Saline(end+1, :) = smoothdata(mean(Lumped_EncoderTrace_R(subjIDX{i}, :), 1), 'gaussian', 240);
    subjMEAN_EncoderValue_P_Saline(end+1, 1) = mean(Lumped_EncoderValue_P(IDX));
    subjMEAN_EncoderValue_R_Saline(end+1, 1) = mean(Lumped_EncoderValue_R(IDX));
    subjMEAN_EncoderTrace_Saline(end+1, :) = mean([Lumped_EncoderTrace_P(IDX, :); Lumped_EncoderTrace_R(IDX, :)], 1);
    subjMEAN_EncoderValue_Saline(end+1, 1) = mean(sessMEAN_EncoderValue_Saline(IDX));
end

sessMEAN_EncoderTrace_P_Saline = Lumped_EncoderTrace_P;
% sessMEAN_EncoderTrace_P_Saline = smoothdata(Lumped_EncoderTrace_P, 2, 'gaussian', 240);
sessMEAN_EncoderTrace_R_Saline = Lumped_EncoderTrace_R;
% sessMEAN_EncoderTrace_R_Saline = smoothdata(Lumped_EncoderTrace_R, 2, 'gaussian', 240);
sessMEAN_EncoderValue_P_Saline = Lumped_EncoderValue_P;
sessMEAN_EncoderValue_R_Saline = Lumped_EncoderValue_R;
% for i = 1:length(Lumped_EncoderValue_P)
%     sessMEAN_EncoderTrace_Saline(end+1, :) = mean([Lumped_EncoderTrace_P(i, :); Lumped_EncoderTrace_R(i, :)], 1);
% end

sessMEAN_EncoderValue_Thm = sessMEAN_EncoderValue_Saline;
%% CNO
% WINst_plt = 2*fs+1;
% WINed_plt = 4*fs;
subjMEAN_spikecnt_CNO = [];
sessMEAN_spikecnt_CNO = nanmean([Lumped_spikecnt_P Lumped_spikecnt_R], 2);
for i = 1:length(subjIDX)
    subjMEAN_spikecnt_CNO(end+1, 1) = nanmean([Lumped_spikecnt_P(subjIDX{i}); Lumped_spikecnt_R(subjIDX{i})]);
end

subjMEAN_EncoderTrace_P_CNO = [];
subjMEAN_EncoderTrace_R_CNO = [];
subjMEAN_EncoderValue_P_CNO = [];
subjMEAN_EncoderValue_R_CNO = [];
subjMEAN_EncoderTrace_CNO = [];
subjMEAN_EncoderValue_CNO = [];
% sessMEAN_EncoderValue_CNO = nanmean([Lumped_EncoderValue_P Lumped_EncoderValue_R], 2);
sessMEAN_EncoderValue_CNO = Lumped_EncoderValue; % comment if only considering [-5s, 0]OC encoder values
for i = 1:length(subjIDX)
    subjMEAN_EncoderTrace_P_CNO(end+1, :) = mean(Lumped_EncoderTrace_P(subjIDX{i}, :), 1);
%     subjMEAN_EncoderTrace_P_CNO(end+1, :) = smoothdata(mean(Lumped_EncoderTrace_P(subjIDX{i}, :), 1), 'gaussian', 240);
    subjMEAN_EncoderTrace_R_CNO(end+1, :) = mean(Lumped_EncoderTrace_R(subjIDX{i}, :), 1);
%     subjMEAN_EncoderTrace_R_CNO(end+1, :) = smoothdata(mean(Lumped_EncoderTrace_R(subjIDX{i}, :), 1), 'gaussian', 240);
    subjMEAN_EncoderValue_P_CNO(end+1, 1) = mean(Lumped_EncoderValue_P(subjIDX{i}));
    subjMEAN_EncoderValue_R_CNO(end+1, 1) = mean(Lumped_EncoderValue_R(subjIDX{i}));
    subjMEAN_EncoderTrace_CNO(end+1, :) = mean([Lumped_EncoderTrace_P(subjIDX{i}, :); Lumped_EncoderTrace_R(subjIDX{i}, :)], 1);
    subjMEAN_EncoderValue_CNO(end+1, 1) = mean(sessMEAN_EncoderValue_CNO(subjIDX{i}));
end

sessMEAN_EncoderTrace_P_CNO = Lumped_EncoderTrace_P;
% sessMEAN_EncoderTrace_P_CNO = smoothdata(Lumped_EncoderTrace_P, 2, 'gaussian', 240);
sessMEAN_EncoderTrace_R_CNO = Lumped_EncoderTrace_R;
% sessMEAN_EncoderTrace_R_CNO = smoothdata(Lumped_EncoderTrace_R, 2, 'gaussian', 240);
sessMEAN_EncoderTrace_CNO = [];
sessMEAN_EncoderValue_P_CNO = Lumped_EncoderValue_P;
sessMEAN_EncoderValue_R_CNO = Lumped_EncoderValue_R;
for i = 1:length(Lumped_EncoderValue_P)
    sessMEAN_EncoderTrace_CNO(end+1, :) = mean([Lumped_EncoderTrace_P(i, :); Lumped_EncoderTrace_R(i, :)], 1);
end
%% Plot average encoder value [-5s 0] by outcome across different brain regions
figure; hold on
bar(1, nanmean(sessMEAN_EncoderValue_OFC))
bar(2, nanmean(sessMEAN_EncoderValue_PPC))
bar(3, nanmean(sessMEAN_EncoderValue_S1))
bar(4, nanmean(sessMEAN_EncoderValue_Thm))
scatter(ones(length(sessMEAN_EncoderValue_OFC), 1), sessMEAN_EncoderValue_OFC)
scatter(2*ones(length(sessMEAN_EncoderValue_PPC), 1), sessMEAN_EncoderValue_PPC)
scatter(3*ones(length(sessMEAN_EncoderValue_S1), 1), sessMEAN_EncoderValue_S1)
scatter(4*ones(length(sessMEAN_EncoderValue_Thm), 1), sessMEAN_EncoderValue_Thm)
errorbar([1 2 3 4], [nanmean(sessMEAN_EncoderValue_OFC) nanmean(sessMEAN_EncoderValue_PPC) nanmean(sessMEAN_EncoderValue_S1) nanmean(sessMEAN_EncoderValue_Thm)], ...
    [SEM(sessMEAN_EncoderValue_OFC) SEM(sessMEAN_EncoderValue_PPC) SEM(sessMEAN_EncoderValue_S1) SEM(sessMEAN_EncoderValue_Thm)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

%% Plot spike rate scatter plot of all sessions
figure; hold on
% subplot(3, 2, 5); hold on
MEAN_spikecnt_P = [];
MEAN_spikecnt_R = [];
for i = 1:length(subjIDX)
    IDX = subjIDX{i}(4:end);
    COLOR = subj_COLOR{i};
    scatter(Lumped_spikecnt_P(IDX), Lumped_spikecnt_R(IDX), 16, 'o', 'MarkerEdgeColor', COLOR)
    MEAN_spikecnt_P(end+1, 1) = nanmean(Lumped_spikecnt_P(IDX));
    MEAN_spikecnt_R(end+1, 1) = nanmean(Lumped_spikecnt_R(IDX));
end
%%%%%%% Add the session index to the dots %%%%%%%
% for i = 1:length(OFCIDX)
%     text(Lumped_NEAch_atOutcome_P(i), Lumped_NEAch_atOutcome_R(i), num2str(i), 'FontSize', 10, 'VerticalAlignment','middle', 'HorizontalAlignment','center')
% end
%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter(nanmean(Lumped_spikecnt_P), nanmean(Lumped_spikecnt_R), 8, 'diamond', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
errorbar(nanmean(Lumped_spikecnt_P), nanmean(Lumped_spikecnt_R), SEM(Lumped_spikecnt_R), 'vertical', 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1)
errorbar(nanmean(Lumped_spikecnt_P), nanmean(Lumped_spikecnt_R), SEM(Lumped_spikecnt_P), 'horizontal', 'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1)
plot([0 nanmean(Lumped_spikecnt_P)], [0 nanmean(Lumped_spikecnt_R)], '-', 'Color', 'r', 'LineWidth' , 0.5)
xlabel('Switch rate before punishment (Hz)')
ylabel('Switch rate before success (Hz)')
Xmin = 0; Xmax = 0.5;
plot([Xmin Xmax], [Xmin Xmax], ':k')
xlim([Xmin Xmax]); ylim([Xmin Xmax])
title([num2str(bandSensor(1)) '~' num2str(bandSensor(2)) 'Hz ' ' ' num2str(WINst_spike) 's~0s before OC     Fs=' num2str(fs) 'Hz'])
text(0.4, 0.3, ['Ratio = ' num2str(nanmean(Lumped_spikecnt_P)/nanmean(Lumped_spikecnt_R))], 'FontSize', 12)
[h, p, t] = qw_statPairedTest(Lumped_spikecnt_P, Lumped_spikecnt_R)
text(0.4, 0.2, ['p < ' num2str(p)], 'FontSize', 12)

figure; hold on
% subplot(3, 2, 6); hold on
scatter(MEAN_spikecnt_P, MEAN_spikecnt_R, 20, 'square', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
Xmin = 0; Xmax = 0.5;
plot([Xmin Xmax], [Xmin Xmax], ':k')
xlim([Xmin Xmax]); ylim([Xmin Xmax])

%% Plot mean spike rate
% Plot mean individual spike rate
figure; hold on
bar(1, mean(subjMEAN_spikecnt_Saline))
bar(2, mean(subjMEAN_spikecnt_CNO))
for i = 1:length(subjMEAN_spikecnt_Saline)
    plot([1 2], [subjMEAN_spikecnt_Saline(i) subjMEAN_spikecnt_CNO(i)], '-ok')
end
errorbar([1 2], [mean(subjMEAN_spikecnt_Saline) mean(subjMEAN_spikecnt_CNO)], [SEM(subjMEAN_spikecnt_Saline) SEM(subjMEAN_spikecnt_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Spike rate (Hz)')
% Plot mean session spike rate
figure; hold on
bar(1, mean(sessMEAN_spikecnt_Saline))
bar(2, mean(sessMEAN_spikecnt_CNO))
scatter(ones(length(sessMEAN_spikecnt_Saline), 1), sessMEAN_spikecnt_Saline)
scatter(2*ones(length(sessMEAN_spikecnt_CNO), 1), sessMEAN_spikecnt_CNO)
errorbar([1 2], [mean(sessMEAN_spikecnt_Saline) mean(sessMEAN_spikecnt_CNO)], [SEM(sessMEAN_spikecnt_Saline) SEM(sessMEAN_spikecnt_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Spike rate (Hz)')
%% Plot spike rate trace before outcomes
winsz = 0.5; step_size = winsz; n = 2; stdev = n/2/4; 
edges1 = -7+winsz/2:step_size:0-winsz/2;
Lumped_FRtrace_P = [];
Lumped_FRtrace_R = [];
for i = 1:len
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(session_spikes_P{i}, step_size, winsz, n, stdev, 10);
    Lumped_FRtrace_P(end+1, :) = fr1(5:14);
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(session_spikes_R{i}, step_size, winsz, n, stdev, 10);
    Lumped_FRtrace_R(end+1, :) = fr1(5:14);
end
figure; hold on
ShadedPlot(edges1, mean(Lumped_FRtrace_P, 1), 'k', 1, SEM(Lumped_FRtrace_P), [0.8 0.8 0.8])
ShadedPlot(edges1, mean(Lumped_FRtrace_R, 1), 'k', 1, SEM(Lumped_FRtrace_R), [0.8 0.8 0.8])
plot(edges1, mean(Lumped_FRtrace_P, 1), '--k', 'LineWidth', 1)
plot(edges1, mean(Lumped_FRtrace_R, 1), '-k', 'LineWidth', 1)
xlabel('Time after Outcome (s)'); ylabel('Switching rate'); vline(0, '-k')
%%
% num_wins = size(Lumped_FRtrace_P, 2)-1;                             
DyFlow_NEAch_atOutcome_P = mean(Lumped_FRtrace_P, 1);
DyFlow_NEAch_atOutcome_R = mean(Lumped_FRtrace_R, 1);
X = DyFlow_NEAch_atOutcome_P(2:end);
Y = DyFlow_NEAch_atOutcome_R(2:end);
figure; hold on
% Create colormap
colormap_name = 'jet';  % Change the colormap name as desired
cmap = colormap(colormap_name);
cmap_length = size(cmap, 1);
% Calculate colors for each vector
num_vectors = size(Lumped_FRtrace_P, 2)-1;
vector_colors = round(linspace(1, cmap_length, num_vectors));
% Plot the vectors with different colors and filled vertices
for i = 2:num_vectors-1 % Plot the vectors with different colors and filled vertices
    x = [X(i), X(i+1)];
    y = [Y(i), Y(i+1)];
    line(x, y, 'Color', cmap(vector_colors(i), :), 'LineWidth', 1.5);
    if i<11
        scatter(x, y, 20, 'filled', 'MarkerEdgeColor', cmap(vector_colors(i), :), 'MarkerFaceColor', cmap(vector_colors(i), :));
        scatter(x, y, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 0.5);
    else
        scatter(x, y, 10, 'filled', 'MarkerEdgeColor', cmap(vector_colors(i), :), 'MarkerFaceColor', cmap(vector_colors(i), :));
    end
end
xlabel('Switching rate before punishment'); ylabel('Switching rate before success');
Xmin = 0; Xmax = 0.4;
plot([Xmin Xmax], [Xmin Xmax], ':k')
xlim([Xmin Xmax]); ylim([Xmin Xmax])
%%%%%% customize colorbar ticks %%%%%%%%
cbh = colorbar;
ylabel(cbh, 'Time after Outcome (s)');
cbh.Limits = [0 1];
% reso = 6;
newTicksLabels = [-5 -4 -3 -2 -1 0 1 2 3];
xst = -5; xed = 3;
slope = (1 - 0)/(xed - xst);
intercept = 0-slope*xst;
newTickValues = newTicksLabels*slope+intercept;
cbh.Ticks = newTickValues;
% a = -10/(num_wins-1); b = 10/(num_wins-1);
cbh.TickLabels = num2cell(newTicksLabels);

% cbh.Direction = 'reverse';



%% Plot mean encoder value (not segregated)
% Plot mean individual encoder value
% figure; hold on
% bar(1, mean(subjMEAN_EncoderValue_Saline))
% % bar(2, mean(subjMEAN_EncoderValue_CNO))
% for i = 1:length(subjMEAN_EncoderValue_Saline)
%     plot([1 2], [subjMEAN_EncoderValue_Saline(i) subjMEAN_EncoderValue_CNO(i)], '-ok')
% end
% errorbar([1 2], [mean(subjMEAN_EncoderValue_Saline) mean(subjMEAN_EncoderValue_CNO)], [SEM(subjMEAN_EncoderValue_Saline) SEM(subjMEAN_EncoderValue_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Mean encoder value (a.u.)')
% Plot mean session spike rate
figure; hold on
bar(1, mean(sessMEAN_EncoderValue_Saline))
% bar(2, mean(sessMEAN_EncoderValue_CNO))
scatter(ones(length(sessMEAN_EncoderValue_Saline), 1), sessMEAN_EncoderValue_Saline)
% scatter(2*ones(length(sessMEAN_EncoderValue_CNO), 1), sessMEAN_EncoderValue_CNO)
% errorbar([1 2], [mean(sessMEAN_EncoderValue_Saline) mean(sessMEAN_EncoderValue_CNO)], [SEM(sessMEAN_EncoderValue_Saline) SEM(sessMEAN_EncoderValue_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
errorbar(1, mean(sessMEAN_EncoderValue_Saline), SEM(sessMEAN_EncoderValue_Saline), 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean encoder value (a.u.)')
%% Plot mean encoder trace (not segregated)
% Plot mean individual encoder trace
TT2 = (WINst_sync*fs:WINed_sync*fs)/120;
figure; hold on
ShadedPlot(TT2, mean(subjMEAN_EncoderTrace_Saline, 1), 'k', 0.5, SEM(subjMEAN_EncoderTrace_Saline), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(subjMEAN_EncoderTrace_CNO, 1), 'k', 0.5, SEM(subjMEAN_EncoderTrace_CNO), [0.9 0.8 0.7])
plot(TT2, mean(subjMEAN_EncoderTrace_Saline, 1), 'LineWidth', 1)
plot(TT2, mean(subjMEAN_EncoderTrace_CNO, 1), 'LineWidth', 1)
vline(0, '-k')
ylabel('Encoder value (a.u.)')
xlabel('Time after outcome (s)')
% Plot mean session encoder trace
figure; hold on
ShadedPlot(TT2, mean(sessMEAN_EncoderTrace_Saline, 1), 'k', 0.5, SEM(sessMEAN_EncoderTrace_Saline), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(sessMEAN_EncoderTrace_CNO, 1), 'k', 0.5, SEM(sessMEAN_EncoderTrace_CNO), [0.9 0.8 0.7])
plot(TT2, mean(sessMEAN_EncoderTrace_Saline, 1), 'LineWidth', 1)
plot(TT2, mean(sessMEAN_EncoderTrace_CNO, 1), 'LineWidth', 1)
vline(0, '-k')
ylabel('Encoder value (a.u.)')
xlabel('Time after outcome (s)')
%% Plot mean encoder value (by each outcome)
% Plot mean individual encoder value

figure; hold on
% subplot(1, 2, 1); hold on
bar(1, mean(subjMEAN_EncoderValue_P_Saline))
bar(2, mean(subjMEAN_EncoderValue_R_Saline))
for i = 1:length(subjMEAN_EncoderValue_P_Saline)
    plot([1 2], [subjMEAN_EncoderValue_P_Saline(i) subjMEAN_EncoderValue_R_Saline(i)], '-ok')
end
errorbar([1 2], [mean(subjMEAN_EncoderValue_P_Saline) mean(subjMEAN_EncoderValue_R_Saline)], [SEM(subjMEAN_EncoderValue_P_Saline) SEM(subjMEAN_EncoderValue_R_Saline)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean encoder value (a.u.)')
% subplot(1, 2, 2); hold on
% bar(1, mean(subjMEAN_EncoderValue_P_CNO))
% bar(2, mean(subjMEAN_EncoderValue_R_CNO))
% for i = 1:length(subjMEAN_EncoderValue_P_CNO)
%     plot([1 2], [subjMEAN_EncoderValue_P_CNO(i) subjMEAN_EncoderValue_R_CNO(i)], '-ok')
% end
% errorbar([1 2], [mean(subjMEAN_EncoderValue_P_CNO) mean(subjMEAN_EncoderValue_R_CNO)], [SEM(subjMEAN_EncoderValue_P_CNO) SEM(subjMEAN_EncoderValue_R_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% Plot mean session encoder value
figure; hold on
% subplot(1, 2, 1); hold on
bar(1, mean(sessMEAN_EncoderValue_P_Saline))
bar(2, mean(sessMEAN_EncoderValue_R_Saline))
scatter(ones(length(sessMEAN_EncoderValue_P_Saline), 1), sessMEAN_EncoderValue_P_Saline)
scatter(2*ones(length(sessMEAN_EncoderValue_R_Saline), 1), sessMEAN_EncoderValue_R_Saline)
errorbar([1 2], [mean(sessMEAN_EncoderValue_P_Saline) mean(sessMEAN_EncoderValue_R_Saline)], [SEM(sessMEAN_EncoderValue_P_Saline) SEM(sessMEAN_EncoderValue_R_Saline)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean encoder value (a.u.)')
% subplot(1, 2, 2); hold on
% bar(1, mean(sessMEAN_EncoderValue_P_CNO))
% bar(2, mean(sessMEAN_EncoderValue_R_CNO))
% scatter(ones(length(sessMEAN_EncoderValue_P_CNO), 1), sessMEAN_EncoderValue_P_CNO)
% scatter(2*ones(length(sessMEAN_EncoderValue_R_CNO), 1), sessMEAN_EncoderValue_R_CNO)
% errorbar([1 2], [mean(sessMEAN_EncoderValue_P_CNO) mean(sessMEAN_EncoderValue_R_CNO)], [SEM(sessMEAN_EncoderValue_P_CNO) SEM(sessMEAN_EncoderValue_R_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
%% Plot mean encoder trace (by each outcome)
% Plot mean individual encoder trace
TT2 = (WINst_sync*fs:WINed_sync*fs)/fs;
figure; hold on
% subplot(1, 2, 1); hold on
% for i = 1:size(subjMEAN_EncoderTrace_P_Saline, 1)
%     plot(TT2, subjMEAN_EncoderTrace_P_Saline(i, :), '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
%     plot(TT2, subjMEAN_EncoderTrace_R_Saline(i, :), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
% end
% ShadedPlot(TT2, mean(subjMEAN_EncoderTrace_P_Saline, 1), 'k', 0.5, SEM(subjMEAN_EncoderTrace_P_Saline), [0.9 0.8 0.7])
% ShadedPlot(TT2, mean(subjMEAN_EncoderTrace_R_Saline, 1), 'k', 0.5, SEM(subjMEAN_EncoderTrace_R_Saline), [0.9 0.8 0.7])
% plot(TT2, mean(subjMEAN_EncoderTrace_P_Saline, 1), '--', 'LineWidth', 1)
% plot(TT2, mean(subjMEAN_EncoderTrace_R_Saline, 1), '-', 'LineWidth', 1)
% ylabel('Encoder value (a.u.)')
% xlabel('Time after outcome (s)')
% vline(0, '-k')
% subplot(1, 2, 2); hold on
% for i = 1:size(subjMEAN_EncoderTrace_P_CNO, 1)
%     plot(TT2, subjMEAN_EncoderTrace_P_CNO(i, :), '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
%     plot(TT2, subjMEAN_EncoderTrace_R_CNO(i, :), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5'

% end
ShadedPlot(TT2, mean(subjMEAN_EncoderTrace_P_CNO, 1), 'k', 0.5, SEM(subjMEAN_EncoderTrace_P_CNO), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(subjMEAN_EncoderTrace_R_CNO, 1), 'k', 0.5, SEM(subjMEAN_EncoderTrace_R_CNO), [0.9 0.8 0.7])
plot(TT2, mean(subjMEAN_EncoderTrace_P_CNO, 1), '--', 'LineWidth', 1)
plot(TT2, mean(subjMEAN_EncoderTrace_R_CNO, 1), '-', 'LineWidth', 1)
xlabel('Time after outcome (s)')
vline(0, '-k')
% Plot mean session encoder trace
figure; hold on
% subplot(1, 2, 1); hold on
% for i = 1:size(sessMEAN_EncoderTrace_P_Saline, 1)
%     plot(TT2, sessMEAN_EncoderTrace_P_Saline(i, :), '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
%     plot(TT2, sessMEAN_EncoderTrace_R_Saline(i, :), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
% end
% ShadedPlot(TT2, mean(sessMEAN_EncoderTrace_P_Saline, 1), 'k', 0.5, SEM(sessMEAN_EncoderTrace_P_Saline), [0.9 0.8 0.7])
% ShadedPlot(TT2, mean(sessMEAN_EncoderTrace_R_Saline, 1), 'k', 0.5, SEM(sessMEAN_EncoderTrace_R_Saline), [0.9 0.8 0.7])
% plot(TT2, mean(sessMEAN_EncoderTrace_P_Saline, 1), '--', 'LineWidth', 1)
% plot(TT2, mean(sessMEAN_EncoderTrace_R_Saline, 1), '-', 'LineWidth', 1)
% ylabel('Encoder value (a.u.)')
% xlabel('Time after outcome (s)')
% vline(0, '-k')
% subplot(1, 2, 2); hold on
% for i = 1:size(sessMEAN_EncoderTrace_P_CNO, 1)
%     plot(TT2, sessMEAN_EncoderTrace_P_CNO(i, :), '--', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
%     plot(TT2, sessMEAN_EncoderTrace_R_CNO(i, :), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5)
% end
ShadedPlot(TT2, mean(sessMEAN_EncoderTrace_P_CNO, 1), 'k', 0.5, SEM(sessMEAN_EncoderTrace_P_CNO), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(sessMEAN_EncoderTrace_R_CNO, 1), 'k', 0.5, SEM(sessMEAN_EncoderTrace_R_CNO), [0.9 0.8 0.7])
plot(TT2, mean(sessMEAN_EncoderTrace_P_CNO, 1), '--', 'LineWidth', 1)
plot(TT2, mean(sessMEAN_EncoderTrace_R_CNO, 1), '-', 'LineWidth', 1)
xlabel('Time after outcome (s)')
vline(0, '-k')

%%
i = 12;
A = session_sync_P{i}(:, 1:601);
B = session_sync_R{i}(:, 1:601);

% Define the time axis
time_axis = (-600:0) / 120;
% Calculate the mean values for the window -3s to -1s
mean_A = mean(A(:, time_axis >= -3 & time_axis <= -1), 2);
mean_B = mean(B(:, time_axis >= -3 & time_axis <= -1), 2);
% Sort the observation rows based on the mean values
[sorted_mean_A, sorted_index_A] = sort(mean_A, 'ascend');
[sorted_mean_B, sorted_index_B] = sort(mean_B, 'ascend');
% Create a figure and subplots
figure;
subplot(2, 1, 1);
imagesc(time_axis, 1:size(A, 1), A(sorted_index_A, :));
colormap(hot);
colorbar;
% vline([-3 -1], '--m')
% xlabel('Time after outcome onset (s)');
title('Failed Trials');

subplot(2, 1, 2);
imagesc(time_axis, 1:size(B, 1), B(sorted_index_B, :));
colormap(hot);
colorbar;
% vline([-3 -1], '--m')
xlabel('Time after outcome onset (s)');
title('Successful Trials');


%% Average session trace of sensor phase, all-trial-lumped sensor phase, session trial-stitched sensor phase by each outcome
WINst_plt = 0*120+1;
WINed_plt = 8*120;
Lumped_NE_phase_P = NaN(len, WINed_plt-WINst_plt+1);
Lumped_Ach_phase_P = NaN(len, WINed_plt-WINst_plt+1);
Lumped_NE_phase_R = NaN(len, WINed_plt-WINst_plt+1);
Lumped_Ach_phase_R = NaN(len, WINed_plt-WINst_plt+1);
for i = 1:len
    Lumped_NE_phase_P(i, :) = nanmean(session_NE_phase_P{i}(:, WINst_plt:WINed_plt), 1);
    Lumped_Ach_phase_P(i, :) = nanmean(session_Ach_phase_P{i}(:, WINst_plt:WINed_plt), 1);
    Lumped_NE_phase_R(i, :) = nanmean(session_NE_phase_R{i}(:, WINst_plt:WINed_plt), 1);
    Lumped_Ach_phase_R(i, :) = nanmean(session_Ach_phase_R{i}(:, WINst_plt:WINed_plt), 1);
end

WINst = 2*120+1;
WINed = 4*120;
LumpedAllTrial_NE_phase_P = [];
LumpedAllTrial_Ach_phase_P = [];
LumpedAllTrial_NE_phase_R = [];
LumpedAllTrial_Ach_phase_R = [];
Lumped_NE_phaseSTITCH_P = cell(len, 1);
Lumped_Ach_phaseSTITCH_P = cell(len, 1);
Lumped_NE_phaseSTITCH_R = cell(len, 1);
Lumped_Ach_phaseSTITCH_R = cell(len, 1);
for i = 1:len
    submatrx = session_NE_phase_P{i}(:, WINst:WINed);
    LumpedAllTrial_NE_phase_P = [LumpedAllTrial_NE_phase_P; submatrx(:)];
    submatrx = session_Ach_phase_P{i}(:, WINst:WINed);
    LumpedAllTrial_Ach_phase_P = [LumpedAllTrial_Ach_phase_P; submatrx(:)];
    submatrx = session_NE_phase_R{i}(:, WINst:WINed);
    LumpedAllTrial_NE_phase_R = [LumpedAllTrial_NE_phase_R; submatrx(:)];
    submatrx = session_Ach_phase_R{i}(:, WINst:WINed);
    LumpedAllTrial_Ach_phase_R = [LumpedAllTrial_Ach_phase_R; submatrx(:)];
    Lumped_NE_phaseSTITCH_P{i} = reshape(session_NE_phase_P{i}(:, WINst:WINed)', 1, []);
    Lumped_Ach_phaseSTITCH_P{i} = reshape(session_Ach_phase_P{i}(:, WINst:WINed)', 1, []);
    Lumped_NE_phaseSTITCH_R{i} = reshape(session_NE_phase_R{i}(:, WINst:WINed)', 1, []);
    Lumped_Ach_phaseSTITCH_R{i} = reshape(session_Ach_phase_R{i}(:, WINst:WINed)', 1, []);
end
%% Plot sensor phase distribution
num_bins = 20;
edges1 = -pi:2*pi/num_bins:pi;
Lumped_NE_phase_distr_P = NaN(len, num_bins);
Lumped_Ach_phase_distr_P = NaN(len, num_bins);
Lumped_NE_phase_distr_R = NaN(len, num_bins);
Lumped_Ach_phase_distr_R = NaN(len, num_bins);
for i = 1:len
    [N, edges] = histcounts(Lumped_NE_phaseSTITCH_P{i}, edges1, 'Normalization', 'probability');
    Lumped_NE_phase_distr_P(i, :) = N;
    [N, edges] = histcounts(Lumped_Ach_phaseSTITCH_P{i}, edges1, 'Normalization', 'probability');
    Lumped_Ach_phase_distr_P(i, :) = N;
    [N, edges] = histcounts(Lumped_NE_phaseSTITCH_R{i}, edges1, 'Normalization', 'probability');
    Lumped_NE_phase_distr_R(i, :) = N;
    [N, edges] = histcounts(Lumped_Ach_phaseSTITCH_R{i}, edges1, 'Normalization', 'probability');
    Lumped_Ach_phase_distr_R(i, :) = N;
end
MEAN_NE_phase_distr_P = [];
MEAN_Ach_phase_distr_P = [];
MEAN_NE_phase_distr_R = [];
MEAN_Ach_phase_distr_R = [];
for i = 1:length(subjIDX)
    MEAN_NE_phase_distr_P(end+1, :) = mean(Lumped_NE_phase_distr_P(subjIDX{i}, :), 1);
    MEAN_Ach_phase_distr_P(end+1, :) = mean(Lumped_Ach_phase_distr_P(subjIDX{i}, :), 1);
    MEAN_NE_phase_distr_R(end+1, :) = mean(Lumped_NE_phase_distr_R(subjIDX{i}, :), 1);
    MEAN_Ach_phase_distr_R(end+1, :) = mean(Lumped_Ach_phase_distr_R(subjIDX{i}, :), 1);
end
figure; hold on
subplot(2, 1, 1); hold on
ShadedPlot(rad2deg(edges(1:end-1)), mean(MEAN_NE_phase_distr_P, 1), 'b', 1, SEM(MEAN_NE_phase_distr_P), [0.9 0.8 0.7])
ShadedPlot(rad2deg(edges(1:end-1)), mean(MEAN_NE_phase_distr_R, 1), 'b', 1, SEM(MEAN_NE_phase_distr_R), [0.9 0.8 0.7])
plot(rad2deg(edges(1:end-1)), mean(MEAN_NE_phase_distr_P, 1), '--b', 'LineWidth', 1)
plot(rad2deg(edges(1:end-1)), mean(MEAN_NE_phase_distr_R, 1), '-b', 'LineWidth', 1)
ylabel('Fraction');
subplot(2, 1, 2); hold on
ShadedPlot(rad2deg(edges(1:end-1)), mean(MEAN_Ach_phase_distr_P, 1), 'r', 1, SEM(MEAN_Ach_phase_distr_P), [0.9 0.8 0.7])
ShadedPlot(rad2deg(edges(1:end-1)), mean(MEAN_Ach_phase_distr_R, 1), 'r', 1, SEM(MEAN_Ach_phase_distr_R), [0.9 0.8 0.7])
plot(rad2deg(edges(1:end-1)), mean(MEAN_Ach_phase_distr_P, 1), '--r', 'LineWidth', 1)
plot(rad2deg(edges(1:end-1)), mean(MEAN_Ach_phase_distr_R, 1), '-r', 'LineWidth', 1)
ylabel('Fraction'); xlabel('Neurotransmitter phase angle (degrees)')

%% Plot average sensor phase trace
MEAN_NE_phase_P = [];
MEAN_Ach_phase_P = [];
MEAN_NE_phase_R = [];
MEAN_Ach_phase_R = [];
for i = 1:length(subjIDX)
    MEAN_NE_phase_P(end+1, :) = mean(Lumped_NE_phase_P(subjIDX{i}, :), 1);
    MEAN_Ach_phase_P(end+1, :) = mean(Lumped_Ach_phase_P(subjIDX{i}, :), 1);
    MEAN_NE_phase_R(end+1, :) = mean(Lumped_NE_phase_R(subjIDX{i}, :), 1);
    MEAN_Ach_phase_R(end+1, :) = mean(Lumped_Ach_phase_R(subjIDX{i}, :), 1);
end
figure; hold on
subplot(2, 1, 1); hold on
TT2 = (WINst_plt:WINed_plt)/120-5;
ShadedPlot(TT2, mean(MEAN_NE_phase_P, 1), 'k', 0.5, SEM(MEAN_NE_phase_P), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(MEAN_NE_phase_R, 1), 'k', 0.5, SEM(MEAN_NE_phase_R), [0.9 0.8 0.7])
plot(TT2, mean(MEAN_NE_phase_P, 1), '--b', 'LineWidth', 1)
plot(TT2, mean(MEAN_NE_phase_R, 1), '-b', 'LineWidth', 1)
ylabel('Phase angle (rads/s)'); vline(0, '-k')

subplot(2, 1, 2); hold on
TT2 = (WINst_plt:WINed_plt)/120-5;
ShadedPlot(TT2, mean(MEAN_Ach_phase_P, 1), 'k', 0.5, SEM(MEAN_Ach_phase_P), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(MEAN_Ach_phase_R, 1), 'k', 0.5, SEM(MEAN_Ach_phase_R), [0.9 0.8 0.7])
plot(TT2, mean(MEAN_Ach_phase_P, 1), '--r', 'LineWidth', 1)
plot(TT2, mean(MEAN_Ach_phase_R, 1), '-r', 'LineWidth', 1)
ylabel('Phase angle (rads/s)'); xlabel('Time after Outcome (s)'); vline(0, '-k')

%% Plot sensor phase distribution
MEAN_NA_phase_P = [];
MEAN_NA_phase_R = [];
num_bins = 120;
for i = 1:len
    MIN = round(-pi, 2); MAX = round(pi, 2);
    norm_NE_phaseSTITCH_P = (Lumped_NE_phaseSTITCH_P{i}-MIN)/(MAX-MIN);
    norm_Ach_phaseSTITCH_P = (Lumped_Ach_phaseSTITCH_P{i}-MIN)/(MAX-MIN);
    [MeanINBin, MeanYINBin, bin_ctrs] = BinnedMean2d(norm_NE_phaseSTITCH_P, norm_Ach_phaseSTITCH_P, [], num_bins);
    MEAN_NA_phase_P(end+1, :, :) = MeanINBin/length(norm_NE_phaseSTITCH_P);
    norm_NE_phaseSTITCH_R = (Lumped_NE_phaseSTITCH_R{i}-MIN)/(MAX-MIN);
    norm_Ach_phaseSTITCH_R = (Lumped_Ach_phaseSTITCH_R{i}-MIN)/(MAX-MIN);
    [MeanINBin, MeanYINBin, bin_ctrs] = BinnedMean2d(norm_NE_phaseSTITCH_R, norm_Ach_phaseSTITCH_R, [], num_bins);
    MEAN_NA_phase_R(end+1, :, :) = MeanINBin/length(norm_NE_phaseSTITCH_R);
end
angle_bins = -180:(180+180)/num_bins:180-(180+180)/num_bins;
A = reshape(mean(MEAN_NA_phase_P, 1), size(mean(MEAN_NA_phase_P, 1), 2), size(mean(MEAN_NA_phase_P, 1), 3));
B = reshape(mean(MEAN_NA_phase_R, 1), size(mean(MEAN_NA_phase_R, 1), 2), size(mean(MEAN_NA_phase_R, 1), 3));
% Create the heatmap
figure;
subplot(2, 1, 1)
imagesc(angle_bins, angle_bins, A);
% Set the colormap and color limits for fractions
colormap(jet);  % You can use other colormaps if desired
caxis([0, max([A(:); B(:)])]);  % Set color limits to [0, max fraction value]
% Add colorbar
colorbar;
% Set axis labels and title
xlabel('Ach phase angle'); ylabel('NE phase angle'); title('Failed')
% Set axis ticks and tick labels
xticks(-180:90:180);
yticks(-180:90:180);
% Flip the y-axis
set(gca, 'YDir', 'reverse');
    
subplot(2, 1, 2)
imagesc(angle_bins, angle_bins, B);
% Set the colormap and color limits for fractions
colormap(jet);  % You can use other colormaps if desired
caxis([0, max([A(:); B(:)])]);  % Set color limits to [0, max fraction value]
% Add colorbar
colorbar;
% Set axis labels and title
xlabel('Ach phase angle'); ylabel('NE phase angle'); title('Successful')
% Set axis ticks and tick labels
xticks(-180:90:180);
yticks(-180:90:180);
% Flip the y-axis
set(gca, 'YDir', 'reverse');

%%
figure; hold on




