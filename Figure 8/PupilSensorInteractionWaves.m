[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);

session_sync = cell(len, 1);
session_spiketime_P = cell(len, 1);
session_spiketime_R = cell(len, 1);
session_OutcomeIDX_P = cell(len, 1);
session_OutcomeIDX_R = cell(len, 1);

winsz = 5;
tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 3; 
bandSensor = [0.4 0.8];
WINspkcnt_st = -5;
WINspkcnt_ed = 1;
EDGE_FR = -3.5:0.5:-1;
num_SHF = 1;
fs = 10;
cnt = 0;
cnt1 = 0;


for i = OFCIDX(:).'
    cnt = cnt+1; 
    S2P_IDX = e1.MetaData.Pupil_in_sensor{i};
    if fs == 120
        z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, [0.1 10], 'bandpass'));
%         z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, 3.5, 'low'));
        z_Ach = zscore(Filter(e1.MetaData.Ach_470{i}, 120, 2, [0.1 10], 'bandpass'));
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
    mx = load(MetaDataX_files{i});
    session_sync{cnt} = sync;

    for k = 2:height(mx.MetaDataX)
        if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>5
            if fs==120
                stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Punish_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Tone_Onset(k));
            else
                stOutcome_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Punish_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Tone_Onset(k));
            end
            idx = intersect(find(spikes==1), stOutcome_sensor-5*fs:stOutcome_sensor);
            labels = [];
            for s1 = 1:length(idx)
                cnt1 = cnt1+1;
                labels(end+1, 1) = cnt1;
            end
            session_spiketime_P{cnt} = [session_spiketime_P{cnt}; [idx labels]];
            session_OutcomeIDX_P{cnt}(end+1, :) = stOutcome_sensor-5*fs:stOutcome_sensor;
        end
        if mx.MetaDataX.(OPPTWIN{option})(k)==1
            if fs==120
                stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Reward_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Tone_Onset(k));
            else
                stOutcome_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Reward_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Tone_Onset(k));
            end
            idx = intersect(find(spikes==1), stOutcome_sensor-5*fs:stOutcome_sensor);
            labels = [];
            for s2 = 1:length(idx)
                cnt1 = cnt1+1;
                labels(end+1, 1) = cnt1;
            end
            session_spiketime_R{cnt} = [session_spiketime_R{cnt}; [idx labels]];
            session_OutcomeIDX_R{cnt}(end+1, :) = stOutcome_sensor-5*fs:stOutcome_sensor;
        end
        
    end
end
toc
%% Lumped waveform observations in different scope (session/subject/all)
session_waveforms_P = cell(len, 1);
session_waveforms_R = cell(len, 1);
session_waveforms_P_plt = cell(len, 1);
session_waveforms_R_plt = cell(len, 1);
wavelength = 2;
wavelength_plt = 6;
cnt = 0;
for i = 1:len
    cnt = cnt+1;
    if isempty(session_spiketime_P{i}) || isempty(session_spiketime_R{i})
        continue;
    end
    session_waveforms_P{cnt} = [extract_spike_sequences(session_sync{i}, session_spiketime_P{i}(:, 1), wavelength, fs) session_spiketime_P{i}(:, 2)];
    session_waveforms_R{cnt} = [extract_spike_sequences(session_sync{i}, session_spiketime_R{i}(:, 1), wavelength, fs) session_spiketime_R{i}(:, 2)];
    session_waveforms_P_plt{cnt} = [extract_spike_sequences(session_sync{i}, session_spiketime_P{i}(:, 1), wavelength_plt, fs) session_spiketime_P{i}(:, 2)];
    session_waveforms_R_plt{cnt} = [extract_spike_sequences(session_sync{i}, session_spiketime_R{i}(:, 1), wavelength_plt, fs) session_spiketime_R{i}(:, 2)];
end

subj_waveforms = cell(length(subjIDX), 1);
for i = 1:length(subjIDX)
    subj_waveforms{i} = cell2mat([session_waveforms_P(subjIDX{i}); session_waveforms_R(subjIDX{i})]);
end

all_waveforms = cell2mat([session_waveforms_P; session_waveforms_R]);

%% Generate the PCA clustering matrix 
K = 2;
D = 1;
PCA_idx = [];
%********** PCA on each subject's waveforms **********%
% for i = 1:length(subjIDX)
%     X = subj_waveforms{i}(:, 1:end-1);
%     [PCcoeff, PCvec] = pca(X);
%     X_pca = X*PCvec(:, 1:D); 
%     [idx, centroids] = kmeans(X_pca, K);
%     PCA_idx = [PCA_idx; [idx subj_waveforms{i}(:, end)]];
% end
%********************************************************
%********** PCA on all sessions' waveforms **********%
X = all_waveforms(:, 1:end-1);
[PCcoeff, PCvec] = pca(X);
X_pca = X*PCvec(:, 1:D); 
[idx, centroids] = kmeans(X_pca, K);
PCA_idx = [idx all_waveforms(:, end)];
 
% scatter_colors = rand(K, 3); 
% figure;
% subplot(2, 1, 1); hold on;
% for k = 1:K
%     scatter(X_pca(idx == k, 1), X_pca(idx == k, 2), [], scatter_colors(k, :), 'filled');
% end
% scatter(centroids(:, 1), centroids(:, 2), 100, scatter_colors, 'x'); 
% xlabel('Principal Component 1');
% ylabel('Principal Component 2');
% subplot(2, 1, 2);
% for k = 1:K
%     avg_waveform = mean(X(idx == k, :));
%     plot(avg_waveform, 'Color', scatter_colors(k, :));
%     hold on;
% end
% xlabel('Time Point');
% ylabel('Amplitude');




%% Update the spiketime matrix by including the PCA clustering indices
session_spiketimePCA_P = cell(len, 1);
session_spiketimePCA_R = cell(len, 1);
session_waveformPCA_P = cell(len, 1);
session_waveformPCA_R = cell(len, 1);
session_Xpca_P = cell(len, 1);
session_Xpca_R = cell(len, 1);
for i = 1:len
    pcidx = find(ismember(PCA_idx(:, 2), session_spiketime_P{i}(:, end))); % locate the trial outcomes in the PCA_idx matrix
    session_spiketimePCA_P{i} = [session_spiketime_P{i} PCA_idx(pcidx, 1)];
    session_waveformPCA_P{i} = [session_waveforms_P{i} PCA_idx(pcidx, 1)];
    session_Xpca_P{i} = [X_pca(pcidx, :) PCA_idx(pcidx, 1)];
    pcidx = find(ismember(PCA_idx(:, 2), session_spiketime_R{i}(:, end))); % locate the trial outcomes in the PCA_idx matrix
    session_spiketimePCA_R{i} = [session_spiketime_R{i} PCA_idx(pcidx, 1)];
    session_waveformPCA_R{i} = [session_waveforms_R{i} PCA_idx(pcidx, 1)];
    session_Xpca_R{i} = [X_pca(pcidx, :) PCA_idx(pcidx, 1)];
end
%%
figure; hold on
for cl = 1:K
    MEAN_Xpca_P = [];
    MEAN_Xpca_R = [];
    MEAN_waveformPCA_P = [];
    MEAN_waveformPCA_R = [];
    for i = 1:len
        clidx_P = find(session_Xpca_P{i}(:, end)==cl);
        clidx_R = find(session_Xpca_R{i}(:, end)==cl);
        if ~isempty(clidx_P)
            MEAN_Xpca_P(end+1, :) = nanmean(session_Xpca_P{i}(clidx_P, 1:end-1), 1);
            MEAN_waveformPCA_P(end+1, :) = nanmean(session_waveformPCA_P{i}(clidx_P, 1:end-2), 1);
        end
        if ~isempty(clidx_R)
            MEAN_Xpca_R(end+1, :) = nanmean(session_Xpca_R{i}(clidx_R, 1:end-1), 1);
            MEAN_waveformPCA_R(end+1, :) = nanmean(session_waveformPCA_R{i}(clidx_R, 1:end-2), 1); 
        end
    end
%     subplot(2, 2, 1); hold on
%     scatter(MEAN_Xpca_P(:, 1), MEAN_Xpca_P(:, 2), 'x', 'MarkerEdgeColor', subj_COLOR{cl})
%     scatter(MEAN_Xpca_R(:, 1), MEAN_Xpca_R(:, 2), 'o', 'MarkerEdgeColor', subj_COLOR{cl})
%     xlabel('PC 1'); ylabel('PC 2');
%     subplot(2, 2, 3); hold on
%     histogram([MEAN_Xpca_P(:, 1); MEAN_Xpca_R(:, 1)], 10, 'FaceColor', subj_COLOR{cl})
% %     xlabel('PC 1'); ylabel('Fraction')
%     subplot(2, 2, 2); hold on
%     histogram([MEAN_Xpca_P(:, 2); MEAN_Xpca_R(:, 2)], 10, 'FaceColor', subj_COLOR{cl})
%     xlabel('PC 2'); ylabel('Fraction')
    subplot(2, 2, 4); hold on
    ShadedPlot((-wavelength*120/2:wavelength*120/2)/120, nanmean(MEAN_waveformPCA_P, 1), subj_COLOR{cl}, 2, SEM(MEAN_waveformPCA_P), [0.8 0.8 0.8])
    ShadedPlot((-wavelength*120/2:wavelength*120/2)/120, nanmean(MEAN_waveformPCA_R, 1), subj_COLOR{cl}, 2, SEM(MEAN_waveformPCA_R), [0.8 0.8 0.8])
    plot((-wavelength*120/2:wavelength*120/2)/120, nanmean(MEAN_waveformPCA_P, 1), '--', 'Color', subj_COLOR{cl})
    plot((-wavelength*120/2:wavelength*120/2)/120, nanmean(MEAN_waveformPCA_R, 1), '-', 'Color', subj_COLOR{cl})
    xlabel('Time after switching (s)'); ylabel('Encoder value')
end

%% Switching-triggered pupil, pupil derivative, HPF and LPF pupil average
Y = session_phiHiPupil;

STA = NaN(len, 6*10+1);
STA_P = NaN(len, 6*10+1);
STA_R = NaN(len, 6*10+1);

SW1_STA_P = NaN(len, 6*10+1);
SW1_STA_R = NaN(len, 6*10+1);

SW2_STA_P = NaN(len, 6*10+1);
SW2_STA_R = NaN(len, 6*10+1);

SHF_STA_P = NaN(len, 6*10+1);
SHF_STA_R = NaN(len, 6*10+1);

cnt=  0;
for i = 1:len
    cnt = cnt+1;
    if isempty(Y{cnt})
        continue;
    elseif isempty(session_spiketime_P{cnt}) || isempty(session_spiketime_R{cnt})
        continue;
    else
        sess_STA = [];
        sess_STA_P = [];
        sess_STA_R = [];
        
        sessSW1_STA_P = [];
        sessSW1_STA_R = [];
        
        sessSW2_STA_P = [];
        sessSW2_STA_R = [];
        
        sessSHF_STA_P = [];
        sessSHF_STA_R = [];
        
        for tt0 = 10:length(session_spiketime{cnt})-5
            sess_STA(end+1, :) = Y{cnt}(session_spiketime{cnt}(tt0)-40:session_spiketime{cnt}(tt0)+20);
        end
        for tt1 = 1:size(session_spiketime_P{cnt}, 1)
            sess_STA_P(end+1, :) = Y{cnt}(session_spiketime_P{cnt}(tt1, 1)-40:session_spiketime_P{cnt}(tt1, 1)+20);
%             if session_spiketime_P{cnt}(tt1, 3)==1
%                 sessSW1_STA_P(end+1, :) = Y{cnt}(session_spiketime_P{cnt}(tt1, 1)-40:session_spiketime_P{cnt}(tt1, 1)+20);
%             else
%                 sessSW2_STA_P(end+1, :) = Y{cnt}(session_spiketime_P{cnt}(tt1, 1)-40:session_spiketime_P{cnt}(tt1, 1)+20);
%             end
        end
        for tt2 = 1:length(session_spiketime_R{cnt})
            sess_STA_R(end+1, :) = Y{cnt}(session_spiketime_R{cnt}(tt2, 1)-40:session_spiketime_R{cnt}(tt2, 1)+20);
%             if session_spiketime_R{cnt}(tt2, 3)==1
%                 sessSW1_STA_R(end+1, :) = Y{cnt}(session_spiketime_R{cnt}(tt2, 1)-40:session_spiketime_R{cnt}(tt2, 1)+20);
%             else
%                 sessSW2_STA_R(end+1, :) = Y{cnt}(session_spiketime_R{cnt}(tt2, 1)-40:session_spiketime_R{cnt}(tt2, 1)+20);
%             end
        end
        for tt3 = 1:length(session_SHFspike_P{cnt})
            sessSHF_STA_P(end+1, :) = Y{cnt}(session_SHFspike_P{cnt}(tt3)-40:session_SHFspike_P{cnt}(tt3)+20);
        end
        for tt4 = 1:length(session_SHFspike_R{cnt})
            sessSHF_STA_R(end+1, :) = Y{cnt}(session_SHFspike_R{cnt}(tt4)-40:session_SHFspike_R{cnt}(tt4)+20);
        end
            
        if ~isempty(sess_STA)
            STA(i, :) = nanmean(sess_STA, 1);
        end
        if ~isempty(sess_STA_P)
            STA_P(i, :) = nanmean(sess_STA_P, 1);
        end
        if ~isempty(sess_STA_R)
            STA_R(i, :) = nanmean(sess_STA_R, 1);
        end
        if ~isempty(sessSW1_STA_P)
            SW1_STA_P(i, :) = nanmean(sessSW1_STA_P, 1);
        end
        if ~isempty(sessSW1_STA_R)
            SW1_STA_R(i, :) = nanmean(sessSW1_STA_R, 1);
        end
        if ~isempty(sessSW2_STA_P)
            SW2_STA_P(i, :) = nanmean(sessSW2_STA_P, 1);
        end
        if ~isempty(sessSW2_STA_R)
            SW2_STA_R(i, :) = nanmean(sessSW2_STA_R, 1);
        end
        if ~isempty(sessSHF_STA_P)
            SHF_STA_P(i, :) = nanmean(sessSHF_STA_P, 1);
        end
        if ~isempty(sessSHF_STA_R)
            SHF_STA_R(i, :) = nanmean(sessSHF_STA_R, 1);
        end
    end
end
%% Visualize K=2 clustering result using the first 2 PCs
figure; hold on
% subplot(1, 3, 1); hold on
ShadedPlot((1:61)/10-4, nanmean(STA, 1), [0 0 0], 2, SEM(STA), [0.8 0.9 0.9])
ShadedPlot((1:61)/10-4, nanmean(STA_P, 1), [0 0 0], 2, SEM(STA_P), [0.8 0.9 0.9])
ShadedPlot((1:61)/10-4, nanmean(STA_R, 1), [0 0 0], 2, SEM(STA_R), [0.8 0.9 0.9])
ShadedPlot((1:61)/10-4, nanmean(SHF_STA_P, 1), [0 0 0], 2, SEM(SHF_STA_P), [0.8 0.9 0.9])
ShadedPlot((1:61)/10-4, nanmean(SHF_STA_R, 1), [0 0 0], 2, SEM(SHF_STA_R), [0.8 0.9 0.9])
plot((1:61)/10-4, nanmean(STA, 1), '-k', 'LineWidth', 1.5)
plot((1:61)/10-4, nanmean(STA_P, 1), '--k', 'LineWidth', 1)
plot((1:61)/10-4, nanmean(STA_R, 1), '-k', 'LineWidth', 1)
plot((1:61)/10-4, nanmean(SHF_STA_P, 1), '--r', 'LineWidth', 1)
plot((1:61)/10-4, nanmean(SHF_STA_R, 1), '-r', 'LineWidth', 1)
title('All Switchings before OC'); vline(0, '-k')
% subplot(2, 3, 4); hold on
% ShadedPlot((1:61)/10-4, nanmean(STA_P, 1)-nanmean(SHF_STA_P, 1), [1 0 0], 2, SEM(STA_P), [0.9 0.8 0.7])
% ShadedPlot((1:61)/10-4, nanmean(STA_R, 1)-nanmean(SHF_STA_R, 1), [1 0 0], 2, SEM(STA_R), [0.9 0.8 0.7])
% plot((1:61)/10-4, nanmean(STA_P, 1)-nanmean(SHF_STA_P, 1), '--', 'Color', [0.64,0.08,0.18], 'LineWidth', 0.5)
% plot((1:61)/10-4, nanmean(STA_R, 1)-nanmean(SHF_STA_R, 1), '-', 'Color', [0.64,0.08,0.18], 'LineWidth', 0.5)
% xlabel('Time after switching (s)'); vline(0, '-k')

% subplot(1, 3, 2); hold on
% ShadedPlot((1:61)/10-4, nanmean(SW1_STA_P, 1), [0 0 0], 2, SEM(SW1_STA_P), [0.8 0.9 0.9])
% ShadedPlot((1:61)/10-4, nanmean(SW1_STA_R, 1), [0 0 0], 2, SEM(SW1_STA_R), [0.8 0.9 0.9])
% plot((1:61)/10-4, nanmean(SW1_STA_P, 1), '--k', 'LineWidth', 0.5)
% plot((1:61)/10-4, nanmean(SW1_STA_R, 1), '-k', 'LineWidth', 0.5)
% title('Short Switchings before OC'); vline(0, '-k')
% subplot(2, 3, 5); hold on
% ShadedPlot((1:61)/10-4, nanmean(SW1_STA_P, 1)-nanmean(SHF_STA_P, 1), [1 0 0], 2, SEM(SW1_STA_P), [0.9 0.8 0.7])
% ShadedPlot((1:61)/10-4, nanmean(SW1_STA_R, 1)-nanmean(SHF_STA_R, 1), [1 0 0], 2, SEM(SW1_STA_R), [0.9 0.8 0.7])
% plot((1:61)/10-4, nanmean(SW1_STA_P, 1)-nanmean(SHF_STA_P, 1), '--', 'Color', [0.64,0.08,0.18], 'LineWidth', 0.5)
% plot((1:61)/10-4, nanmean(SW1_STA_R, 1)-nanmean(SHF_STA_R, 1), '-', 'Color', [0.64,0.08,0.18], 'LineWidth', 0.5)
% xlabel('Time after switching (s)'); vline(0, '-k')

% subplot(1, 3, 3); hold on
% ShadedPlot((1:61)/10-4, nanmean(SW2_STA_P, 1), [0 0 0], 2, SEM(SW2_STA_P), [0.8 0.9 0.9])
% ShadedPlot((1:61)/10-4, nanmean(SW2_STA_R, 1), [0 0 0], 2, SEM(SW2_STA_R), [0.8 0.9 0.9])
% plot((1:61)/10-4, nanmean(SW2_STA_P, 1), '--k', 'LineWidth', 2)
% plot((1:61)/10-4, nanmean(SW2_STA_R, 1), '-k', 'LineWidth', 2)
% title('Long Switchings before OC'); vline(0, '-k')
% subplot(2, 3, 6); hold on
% ShadedPlot((1:61)/10-4, nanmean(SW2_STA_P, 1)-nanmean(SHF_STA_P, 1), [1 0 0], 2, SEM(SW2_STA_P), [0.9 0.8 0.7])
% ShadedPlot((1:61)/10-4, nanmean(SW2_STA_R, 1)-nanmean(SHF_STA_R, 1), [1 0 0], 2, SEM(SW2_STA_R), [0.9 0.8 0.7])
% plot((1:61)/10-4, nanmean(SW2_STA_P, 1)-nanmean(SHF_STA_P, 1), '--', 'Color', [0.64,0.08,0.18], 'LineWidth', 0.5)
% plot((1:61)/10-4, nanmean(SW2_STA_R, 1)-nanmean(SHF_STA_R, 1), '-', 'Color', [0.64,0.08,0.18], 'LineWidth', 0.5)
% xlabel('Time after switching (s)'); vline(0, '-k')

%% Spike event series labeled by PC-spike value (1 or 2)
spike_atOutcome_P = cell(len, 1);
spike_atOutcome_R = cell(len, 1);
for i = 1:len
    for j = 1:size(session_OutcomeIDX_P{i}, 1)
        spikeidx_time = find(ismember(session_OutcomeIDX_P{i}(j, :), session_spiketimePCA_P{i}(:, 1)));
        spikeidx_spike = find(ismember(session_spiketimePCA_P{i}(:, 1), session_OutcomeIDX_P{i}(j, :)));
        trl_spikes = zeros(1, length(session_OutcomeIDX_P{i}(j, :)));
        trl_spikes(spikeidx_time) = session_spiketimePCA_P{i}(spikeidx_spike, 3)';
        spike_atOutcome_P{i}(end+1, :) = trl_spikes;
    end
    for j = 1:size(session_OutcomeIDX_R{i}, 1)
        spikeidx_time = find(ismember(session_OutcomeIDX_R{i}(j, :), session_spiketimePCA_R{i}(:, 1)));
        spikeidx_spike = find(ismember(session_spiketimePCA_R{i}(:, 1), session_OutcomeIDX_R{i}(j, :)));
        trl_spikes = zeros(1, length(session_OutcomeIDX_R{i}(j, :)));
        trl_spikes(spikeidx_time) = session_spiketimePCA_R{i}(spikeidx_spike, 3)';
        spike_atOutcome_R{i}(end+1, :) = trl_spikes;
    end
end
%% Calculate the PC-spike rate
spike1_atOutcome_P = spike_atOutcome_P;
spike1_atOutcome_R = spike_atOutcome_R;
spike2_atOutcome_P = spike_atOutcome_P;
spike2_atOutcome_R = spike_atOutcome_R;
for i = 1:len
    for j = 1:size(spike_atOutcome_P{i}, 1)
        spike1_atOutcome_P{i}(j, find(spike_atOutcome_P{i}(j, :)~=1)) = 0;
        spike2_atOutcome_P{i}(j, find(spike_atOutcome_P{i}(j, :)~=2)) = 0;
        spike2_atOutcome_P{i}(j, find(spike_atOutcome_P{i}(j, :)==2)) = 1;
    end
    for j = 1:size(spike_atOutcome_R{i}, 1)
        spike1_atOutcome_R{i}(j, find(spike_atOutcome_R{i}(j, :)~=1)) = 0;
        spike2_atOutcome_R{i}(j, find(spike_atOutcome_R{i}(j, :)~=2)) = 0;
        spike2_atOutcome_R{i}(j, find(spike_atOutcome_R{i}(j, :)==2)) = 1;
    end
end
%%
winsz = 0.5; step_size = winsz; n = 2; stdev = n/2/4; 
edges1 = -5+winsz/2:step_size:0-winsz/2;
Lumped_spike1FR_P = [];
Lumped_spike1FR_R = [];
Lumped_spike2FR_P = [];
Lumped_spike2FR_R = [];
for i = 1:len
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(spike1_atOutcome_P{i}, step_size, winsz, n, stdev, 120);
    Lumped_spike1FR_P(end+1, :) = fr0;
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(spike1_atOutcome_R{i}, step_size, winsz, n, stdev, 120);
    Lumped_spike1FR_R(end+1, :) = fr0;
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(spike2_atOutcome_P{i}, step_size, winsz, n, stdev, 120);
    Lumped_spike2FR_P(end+1, :) = fr0;
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(spike2_atOutcome_R{i}, step_size, winsz, n, stdev, 120);
    Lumped_spike2FR_R(end+1, :) = fr0;
end
%%
figure; hold on
subplot(1, 2, 1); hold on
ShadedPlot(edges1, nanmean(Lumped_spike1FR_P, 1), [0 0 0], 2, SEM(Lumped_spike1FR_P), [0.8 0.8 0.8])
ShadedPlot(edges1, nanmean(Lumped_spike1FR_R, 1), [0 0 0], 2, SEM(Lumped_spike1FR_R), [0.8 0.8 0.8])
plot(edges1, nanmean(Lumped_spike1FR_P, 1), '--', 'Color', subj_COLOR{1}, 'LineWidth', 1)
plot(edges1, nanmean(Lumped_spike1FR_R, 1), '-', 'Color', subj_COLOR{1}, 'LineWidth', 1)
subplot(1, 2, 2); hold on
ShadedPlot(edges1, nanmean(Lumped_spike2FR_P, 1), [0 0 0], 2, SEM(Lumped_spike2FR_P), [0.8 0.8 0.8])
ShadedPlot(edges1, nanmean(Lumped_spike2FR_R, 1), [0 0 0], 2, SEM(Lumped_spike2FR_R), [0.8 0.8 0.8])
plot(edges1, nanmean(Lumped_spike2FR_P, 1), '--', 'Color', subj_COLOR{2}, 'LineWidth', 1)
plot(edges1, nanmean(Lumped_spike2FR_R, 1), '-', 'Color', subj_COLOR{2}, 'LineWidth', 1)



%% Sort waveforms based on the PC1
D = 1;
%********** PCA on each subject's waveforms **********%
% for i = 1:length(subjIDX)
%     X = subj_waveforms{i}(:, 1:end-1);
%     [PCcoeff, PCvec] = pca(X);
%     X_pca = X*PCvec(:, 1:D); 
%     [idx, centroids] = kmeans(X_pca, K);
%     PCA_idx = [PCA_idx; [idx subj_waveforms{i}(:, end)]];
% end
%********************************************************
%********** PCA on all sessions' waveforms **********%
X = all_waveforms(:, 1:end-1);
[PCcoeff, PCvec] = pca(X);
X_pca = X*PCvec(:, 1:D); 
% store waveform (its spike global index) of each quantile
num_qt = 3;
X_pca_idx = [X_pca all_waveforms(:, end)];
qt = 1/num_qt:1/num_qt:1-1/num_qt;
X_PC1qt = cell(1, num_qt);
PC1qt = quantile(X_pca(:, 1),  qt);
PC1qt = [-Inf PC1qt Inf];
for i = 1:num_qt
    X_PC1qt{i} = X_pca_idx(find(X_pca_idx(:, 1)>=PC1qt(i)&X_pca_idx(:, 1)<PC1qt(i+1)), 2);
end
% get session outcome waveforms with its quantile assignment
session_waveformPC1_P = cell(len, 1);
session_waveformPC1_R = cell(len, 1);
for i = 1:len
    label_waveforms_P = NaN(size(session_waveforms_P{i}, 1), 1);
    label_waveforms_R = NaN(size(session_waveforms_R{i}, 1), 1);
    if isempty(session_waveforms_P{i}) || isempty(session_waveforms_R{i})
        continue
    end
    for cl = 1:num_qt
        label_waveforms_P(find(ismember(session_waveforms_P{i}(:, end), X_PC1qt{cl}))) = cl;
        label_waveforms_R(find(ismember(session_waveforms_R{i}(:, end), X_PC1qt{cl}))) = cl;
    end
    session_waveformPC1_P{i} = [UnitNormalization(session_waveforms_P_plt{i}(:, 1:end-1)) label_waveforms_P];
    session_waveformPC1_R{i} = [UnitNormalization(session_waveforms_R_plt{i}(:, 1:end-1)) label_waveforms_R];
end

figure; hold on
for cl = [1 2 3]
    MEAN_waveformPC1_P = [];
    MEAN_waveformPC1_R = [];
    for i = 1:len
        if isempty(session_waveformPC1_P{i}) || isempty(session_waveformPC1_R{i})
            continue
        end
        MEAN_waveformPC1_P(end+1, :) = nanmean(session_waveformPC1_P{i}(find(session_waveformPC1_P{i}(:, end)==cl), 1:end-1), 1);
        MEAN_waveformPC1_R(end+1, :) = nanmean(session_waveformPC1_R{i}(find(session_waveformPC1_R{i}(:, end)==cl), 1:end-1), 1);
    end      
    ShadedPlot((-wavelength_plt*fs/2:wavelength_plt*fs/2)/fs, nanmean(MEAN_waveformPC1_P, 1), subj_COLOR{cl}, 2, SEM(MEAN_waveformPC1_P), [0.8 0.8 0.8])
    ShadedPlot((-wavelength_plt*fs/2:wavelength_plt*fs/2)/fs, nanmean(MEAN_waveformPC1_R, 1), subj_COLOR{cl}, 2, SEM(MEAN_waveformPC1_R), [0.8 0.8 0.8])
    plot((-wavelength_plt*fs/2:wavelength_plt*fs/2)/fs, nanmean(MEAN_waveformPC1_P, 1), '--', 'Color', subj_COLOR{cl})
    plot((-wavelength_plt*fs/2:wavelength_plt*fs/2)/fs, nanmean(MEAN_waveformPC1_R, 1), '-', 'Color', subj_COLOR{cl})
end
xlabel('Time after switching (s)'); ylabel('Encoder value')
%% Calculate the event rate of quantile-sorted waveforms
PC1spike_atOutcome_P = cell(len, 1);
PC1spike_atOutcome_R = cell(len, 1);
for i = 1:len
    if isempty(session_waveformPC1_P{i}) || isempty(session_waveformPC1_R{i})
        continue
    end
    for j = 1:size(session_OutcomeIDX_P{i}, 1)
        spikeidx_time = find(ismember(session_OutcomeIDX_P{i}(j, :), session_spiketime_P{i}(:, 1)));
        spikeidx_spike = find(ismember(session_spiketime_P{i}(:, 1), session_OutcomeIDX_P{i}(j, :)));
        trl_spikes = zeros(1, length(session_OutcomeIDX_P{i}(j, :)));
        trl_spikes(spikeidx_time) = session_waveformPC1_P{i}(spikeidx_spike, end)';
        PC1spike_atOutcome_P{i}(end+1, :) = trl_spikes;
    end
    for j = 1:size(session_OutcomeIDX_R{i}, 1)
        spikeidx_time = find(ismember(session_OutcomeIDX_R{i}(j, :), session_spiketime_R{i}(:, 1)));
        spikeidx_spike = find(ismember(session_spiketime_R{i}(:, 1), session_OutcomeIDX_R{i}(j, :)));
        trl_spikes = zeros(1, length(session_OutcomeIDX_R{i}(j, :)));
        trl_spikes(spikeidx_time) = session_waveformPC1_R{i}(spikeidx_spike, end)';
        PC1spike_atOutcome_R{i}(end+1, :) = trl_spikes;
    end
end

winsz = 0.5; step_size = winsz; n = 2; stdev = n/2/4; 
edges1 = -5+winsz/2:step_size:0-winsz/2;
figure; hold on
for cl = [1 2 3]
    QTspike_atOutcome_P = PC1spike_atOutcome_P;
    QTspike_atOutcome_R = PC1spike_atOutcome_R;
    
    Lumped_PC1spikeFR_P = [];
    Lumped_PC1spikeFR_R = [];
    for i = 1:len
        if isempty(PC1spike_atOutcome_P{i}) || isempty(PC1spike_atOutcome_R{i})
            continue
        end
        for j = 1:size(PC1spike_atOutcome_P{i}, 1)
            QTspike_atOutcome_P{i}(j, find(PC1spike_atOutcome_P{i}(j, :)~=cl)) = 0;
            QTspike_atOutcome_P{i}(j, find(PC1spike_atOutcome_P{i}(j, :)==cl)) = 1;
        end
        for j = 1:size(PC1spike_atOutcome_R{i}, 1)
            QTspike_atOutcome_R{i}(j, find(PC1spike_atOutcome_R{i}(j, :)~=cl)) = 0;
            QTspike_atOutcome_R{i}(j, find(PC1spike_atOutcome_R{i}(j, :)==cl)) = 1;
        end
        
        [fr0, fr1, sdf0, sdf1] = sdfEstimate(QTspike_atOutcome_P{i}, step_size, winsz, n, stdev, fs);
        Lumped_PC1spikeFR_P(end+1, :) = fr0;
        [fr0, fr1, sdf0, sdf1] = sdfEstimate(QTspike_atOutcome_R{i}, step_size, winsz, n, stdev, fs);
        Lumped_PC1spikeFR_R(end+1, :) = fr0;
    end
    subplot(2, num_qt, cl); hold on; title([num2str(cl) '/' num2str(num_qt)])
    ShadedPlot(edges1, nanmean(Lumped_PC1spikeFR_P, 1), [0 0 0], 2, SEM(Lumped_PC1spikeFR_P), [0.8 0.8 0.8])
    ShadedPlot(edges1, nanmean(Lumped_PC1spikeFR_R, 1), [0 0 0], 2, SEM(Lumped_PC1spikeFR_R), [0.8 0.8 0.8])
    plot(edges1, nanmean(Lumped_PC1spikeFR_P, 1), '--', 'Color', subj_COLOR{cl}, 'LineWidth', 1)
    plot(edges1, nanmean(Lumped_PC1spikeFR_R, 1), '-', 'Color', subj_COLOR{cl}, 'LineWidth', 1)
    ylim([0 0.4])
    
    subplot(2, num_qt, 3+cl); hold on
    bar(1, nanmean(nanmean(Lumped_PC1spikeFR_P(:, 5:end), 2), 1))
    bar(2, nanmean(nanmean(Lumped_PC1spikeFR_R(:, 5:end), 2), 1))
    for k = 1:size(Lumped_PC1spikeFR_P, 1)
        plot([1 2], [nanmean(Lumped_PC1spikeFR_P(k, 5:end), 2) nanmean(Lumped_PC1spikeFR_R(k, 5:end), 2)], '-ok')
    end
    errorbar([1 2], [nanmean(nanmean(Lumped_PC1spikeFR_P(:, 5:end), 2), 1) nanmean(nanmean(Lumped_PC1spikeFR_R(:, 5:end), 2), 1)], ...
        [SEM(nanmean(Lumped_PC1spikeFR_P(:, 5:end), 2)) SEM(nanmean(Lumped_PC1spikeFR_R(:, 5:end), 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
    [h, p, t] = qw_statPairedTest(nanmean(Lumped_PC1spikeFR_P(:, 5:end), 2), nanmean(Lumped_PC1spikeFR_R(:, 5:end), 2));
    text(1.5, 0.3, ['p < ' num2str(p)], 'FontSize', 8)
    ylim([0 0.4])
end

%% Calculate the dilation rate 
dila_atOutcome_P = cell(len, 1);
dila_atOutcome_R = cell(len, 1);
for i = 1:len
    if isempty(session_dilaPupil{i})
        continue
    end
    for j = 1:size(session_OutcomeIDX_P{i}, 1)
        dilaidx_time = find(ismember(session_OutcomeIDX_P{i}(j, :), session_dilaPupil{i}));
        trl_spikes = zeros(1, length(session_OutcomeIDX_P{i}(j, :)));
        trl_spikes(dilaidx_time) = 1;
        dila_atOutcome_P{i}(end+1, :) = trl_spikes;
    end
    for j = 1:size(session_OutcomeIDX_R{i}, 1)
        dilaidx_time = find(ismember(session_OutcomeIDX_R{i}(j, :), session_dilaPupil{i}));
        trl_spikes = zeros(1, length(session_OutcomeIDX_R{i}(j, :)));
        trl_spikes(dilaidx_time) = 1;
        dila_atOutcome_R{i}(end+1, :) = trl_spikes;
    end
end

winsz = 0.5; step_size = winsz; n = 2; stdev = n/2/4; 
edges1 = -7+winsz/2:step_size:0-winsz/2;
edges2 = -5+winsz/2:step_size:0-winsz/2;
figure; hold on
    
Lumped_spikeFR_P = [];
Lumped_spikeFR_R = [];
Lumped_dilaFR_P = [];
Lumped_dilaFR_R = [];
for i = 1:len
    if isempty(dila_atOutcome_P{i}) || isempty(dila_atOutcome_R{i})
        continue
    end
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(session_spikes_P{i}, step_size, winsz, n, stdev, fs);
    Lumped_spikeFR_P(end+1, :) = fr1(5:14);
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(session_spikes_R{i}, step_size, winsz, n, stdev, fs);
    Lumped_spikeFR_R(end+1, :) = fr1(5:14);

    [fr0, fr1, sdf0, sdf1] = sdfEstimate(dila_atOutcome_P{i}, step_size, winsz, n, stdev, fs);
    Lumped_dilaFR_P(end+1, :) = fr0;
    [fr0, fr1, sdf0, sdf1] = sdfEstimate(dila_atOutcome_R{i}, step_size, winsz, n, stdev, fs);
    Lumped_dilaFR_R(end+1, :) = fr0;
end

subplot(1, 2, 1); hold on
ShadedPlot(edges2, nanmean(Lumped_spikeFR_P, 1), [0 0 0], 2, SEM(Lumped_spikeFR_P), [0.8 0.8 0.8])
ShadedPlot(edges2, nanmean(Lumped_dilaFR_P, 1), [0 0 0], 2, SEM(Lumped_dilaFR_P), [0.8 0.8 0.8])
plot(edges2, nanmean(Lumped_spikeFR_P, 1), '--', 'Color', [0 0 0], 'LineWidth', 1)
plot(edges2, nanmean(Lumped_dilaFR_P, 1), '--', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1)
subplot(1, 2, 2); hold on
ShadedPlot(edges2, nanmean(Lumped_spikeFR_R, 1), [0 0 0], 2, SEM(Lumped_spikeFR_R), [0.8 0.8 0.8])
ShadedPlot(edges2, nanmean(Lumped_dilaFR_R, 1), [0 0 0], 2, SEM(Lumped_dilaFR_R), [0.8 0.8 0.8])
plot(edges2, nanmean(Lumped_spikeFR_R, 1), '-', 'Color', [0 0 0], 'LineWidth', 1)
plot(edges2, nanmean(Lumped_dilaFR_R, 1), '-', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1)
%%
MEAN_spikexdila_PCC_P = [];
MEAN_spikexdila_PCC_R = [];
for i = 1:size(Lumped_spikeFR_P, 1)
    R = corrcoef(Lumped_spikeFR_P(i, :), Lumped_dilaFR_P(i, :));
    MEAN_spikexdila_PCC_P(end+1, 1) = R(2);
    R = corrcoef(Lumped_spikeFR_R(i, :), Lumped_dilaFR_R(i, :));
    MEAN_spikexdila_PCC_R(end+1, 1) = R(2);
end

figure; hold on
bar(1, nanmean(MEAN_spikexdila_PCC_P))
bar(2, nanmean(MEAN_spikexdila_PCC_R))
for i = 1:length(MEAN_spikexdila_PCC_P)
    plot(1:2, [MEAN_spikexdila_PCC_P(i) MEAN_spikexdila_PCC_R(i)], '-k')
end
errorbar(1:2, [nanmean(MEAN_spikexdila_PCC_P) nanmean(MEAN_spikexdila_PCC_R)], [SEM(MEAN_spikexdila_PCC_P) SEM(MEAN_spikexdila_PCC_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
    
    









