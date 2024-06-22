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

fs = 120;
Lumped_NE_atTone_P = NaN(len, 5*fs+1*fs+1); % Punished but withheld >4s NE
Lumped_Ach_atTone_P = NaN(len, 5*fs+1*fs+1); % Punished but withheld >4s Ach
Lumped_NE_atTone_R = NaN(len, 5*fs+1*fs+1); % Rewarded NE
Lumped_Ach_atTone_R = NaN(len, 5*fs+1*fs+1); % Rewarded Ach
Lumped_NE_atOutcome_P = NaN(len, 5*fs+1*fs+1); % Punished but withheld >4s NE
Lumped_Ach_atOutcome_P = NaN(len, 5*fs+1*fs+1); % Punished but withheld >4s Ach
Lumped_NE_atOutcome_R = NaN(len, 5*fs+1*fs+1); % Rewarded NE
Lumped_Ach_atOutcome_R = NaN(len, 5*fs+1*fs+1); % Rewarded Ach

session_NE_P = cell(len, 1);
session_NE_R = cell(len, 1);
session_Ach_P = cell(len, 1);
session_Ach_R = cell(len, 1);
tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 3;
bandSensor = [0.4 0.8];
HPSensor = 0.1;
fs = 120;
WIN_spikecnt = 4;

cnt = 0;
for i = OFCIDX(:).'
    cnt = cnt+1;
    S2P_IDX = e1.MetaData.Pupil_in_sensor{i};
    if fs == 120
        z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, [HPSensor 5], 'bandpass'));
%         z_NE = zscore(e1.MetaData.NE_470{i});
        z_Ach = zscore(Filter(e1.MetaData.Ach_470{i}, 120, 2, [HPSensor 5], 'bandpass'));
%         z_Ach = zscore(e1.MetaData.Ach_470{i});
        i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass');
        i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass');
        timestamp = e1.MetaData.Sensor_time{i};
    else
        z_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, [HPSensor 4], 'bandpass'); z_NE = zscore(z_NE(S2P_IDX));
%         z_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, 3.5, 'low'); z_NE = zscore(z_NE(S2P_IDX));
        z_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, [HPSensor 4], 'bandpass'); z_Ach = zscore(z_Ach(S2P_IDX));
%         z_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, 3.5, 'low'); z_Ach = zscore(z_Ach(S2P_IDX));
        i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass'); i_NE = zscore(i_NE(S2P_IDX));
        i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass'); i_Ach = zscore(i_Ach(S2P_IDX));
        timestamp = e1.MetaData.Pupil_time{i};
    end

    [phi_NE, phi_Ach, spikes, sync, spikes_01, sync_hilbert] = CalcSpikes(i_NE, i_Ach, fs);
%     m = load(MetaData_files{i});
    mx = load(MetaDataX_files{i});
    session_NE_atTone_P = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_Ach_atTone_P = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_NE_atTone_R = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_Ach_atTone_R = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_NE_atOutcome_P = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_Ach_atOutcome_P = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_NE_atOutcome_R = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    session_Ach_atOutcome_R = NaN(height(mx.MetaDataX), 5*fs+1*fs+1);
    for k = 1:height(mx.MetaDataX)
        if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>5
            stOutcome_sensor = CrossSampling(timestamp, mx.MetaDataX.Punish_Onset(k));
            stTone_sensor = CrossSampling(timestamp, mx.MetaDataX.Tone_Onset(k));
            if stTone_sensor+4*fs <= length(z_NE)
                session_NE_atTone_P(k, :) = z_NE(stTone_sensor-2*fs:stTone_sensor+4*fs);
                session_Ach_atTone_P(k, :) = z_Ach(stTone_sensor-2*fs:stTone_sensor+4*fs);
            end
            if stOutcome_sensor+1*fs <= length(z_NE)
                session_NE_atOutcome_P(k, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
                session_Ach_atOutcome_P(k, :) = z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            end
        end
        if mx.MetaDataX.(OPPTWIN{option})(k)==1
        %%%%%%%%% for new S1 animals %%%%%%%%%%%%
%         if mx.MetaDataX.(OPPTWIN{option})(k)~=0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            stOutcome_sensor = CrossSampling(timestamp, mx.MetaDataX.Reward_Onset(k));
            stTone_sensor = CrossSampling(timestamp, mx.MetaDataX.Tone_Onset(k));
            if stTone_sensor+4*fs <= length(z_NE)
                session_NE_atTone_R(k, :) = z_NE(stTone_sensor-2*fs:stTone_sensor+4*fs);
                session_Ach_atTone_R(k, :) = z_Ach(stTone_sensor-2*fs:stTone_sensor+4*fs);
            end
            if stOutcome_sensor+1*fs <= length(z_NE)
                session_NE_atOutcome_R(k, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
                session_Ach_atOutcome_R(k, :) = z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            end
        end
    end
    Lumped_NE_atTone_P(cnt, :) = nanmean(session_NE_atTone_P, 1); 
    Lumped_Ach_atTone_P(cnt, :) = nanmean(session_Ach_atTone_P, 1); 
    Lumped_NE_atTone_R(cnt, :) = nanmean(session_NE_atTone_R, 1); 
    Lumped_Ach_atTone_R(cnt, :) = nanmean(session_Ach_atTone_R, 1); 
    Lumped_NE_atOutcome_P(cnt, :) = nanmean(session_NE_atOutcome_P, 1); 
    Lumped_Ach_atOutcome_P(cnt, :) = nanmean(session_Ach_atOutcome_P, 1); 
    Lumped_NE_atOutcome_R(cnt, :) = nanmean(session_NE_atOutcome_R, 1); 
    Lumped_Ach_atOutcome_R(cnt, :) = nanmean(session_Ach_atOutcome_R, 1);
    
    session_NE_atOutcome_P(find(isnan(session_NE_atOutcome_P(:, 1))), :) = []; session_NE_P{cnt} = session_NE_atOutcome_P;
    session_NE_atOutcome_R(find(isnan(session_NE_atOutcome_R(:, 1))), :) = []; session_NE_R{cnt} = session_NE_atOutcome_R;
    session_Ach_atOutcome_P(find(isnan(session_Ach_atOutcome_P(:, 1))), :) = []; session_Ach_P{cnt} = session_Ach_atOutcome_P;
    session_Ach_atOutcome_R(find(isnan(session_Ach_atOutcome_R(:, 1))), :) = []; session_Ach_R{cnt} = session_Ach_atOutcome_R;
end
toc

%%
TT1 = (-2*fs:4*fs)/120;
subj_NE_atTone_P = NaN(length(subjIDX), 721);
subj_Ach_atTone_P = NaN(length(subjIDX), 721);
    subj_NE_atTone_R = NaN(length(subjIDX), 721);
    subj_Ach_atTone_R = NaN(length(subjIDX), 721);
TT2 = (-5*fs:1*fs)/fs;
subj_NE_atOutcome_P = NaN(length(subjIDX), length(TT2));
subj_Ach_atOutcome_P = NaN(length(subjIDX), length(TT2));
    subj_NE_atOutcome_R = NaN(length(subjIDX), length(TT2));
    subj_Ach_atOutcome_R = NaN(length(subjIDX), length(TT2));
    
for i = 1:length(subjIDX)
    subj_NE_atTone_P(i, :) = nanmean(Lumped_NE_atTone_P(subjIDX{i}, :), 1);
    subj_Ach_atTone_P(i, :) = nanmean(Lumped_Ach_atTone_P(subjIDX{i}, :), 1);
        subj_NE_atTone_R(i, :) = nanmean(Lumped_NE_atTone_R(subjIDX{i}, :), 1);
        subj_Ach_atTone_R(i, :) = nanmean(Lumped_Ach_atTone_R(subjIDX{i}, :), 1);
    subj_NE_atOutcome_P(i, :) = nanmean(Lumped_NE_atOutcome_P(subjIDX{i}, :), 1);
    subj_Ach_atOutcome_P(i, :) = nanmean(Lumped_Ach_atOutcome_P(subjIDX{i}, :), 1);
        subj_NE_atOutcome_R(i, :) = nanmean(Lumped_NE_atOutcome_R(subjIDX{i}, :), 1);
        subj_Ach_atOutcome_R(i, :) = nanmean(Lumped_Ach_atOutcome_R(subjIDX{i}, :), 1);
end
%%
figure; hold on
subplot(2, 1, 1); hold on
ShadedPlot(TT1, mean(subj_NE_atTone_P, 1), [0 0 1], 1, SEM(subj_NE_atTone_P), [0.73 0.83 0.96])
plot(TT1, mean(subj_NE_atTone_P, 1), '--b', 'LineWidth', 1)
ShadedPlot(TT1, mean(subj_NE_atTone_R, 1), [0 0 1], 1, SEM(subj_NE_atTone_R), [0.73 0.83 0.96])
plot(TT1, mean(subj_NE_atTone_R, 1), '-b', 'LineWidth', 1)
ylabel('Z-score')
vline(0, '-k')
subplot(2, 1, 2); hold on
ShadedPlot(TT1, mean(subj_Ach_atTone_P, 1), [1 0 0], 1, SEM(subj_Ach_atTone_P), [0.9 0.8 0.7])
plot(TT1, mean(subj_Ach_atTone_P, 1), '--r', 'LineWidth', 1)
ShadedPlot(TT1, mean(subj_Ach_atTone_R, 1), [1 0 0], 1, SEM(subj_Ach_atTone_R), [0.9 0.8 0.7])
plot(TT1, mean(subj_Ach_atTone_R, 1), '-r', 'LineWidth', 1)
ylabel('Z-score')
xlabel('Time after Inhibition Tone Onset (s)')
vline(0, '-k')
%%
TT2 = (-5*fs:1*fs)/fs;
%%%% for S1 animals, discard the traces with dips %%%%%%%%%
% for i = 1:size(Lumped_Ach_atOutcome_R, 1)
%     if ~isempty(find(Lumped_Ach_atOutcome_R(i, :)<-2))
%         Lumped_Ach_atOutcome_R(i, :) = NaN(1, size(Lumped_Ach_atOutcome_R, 2));
%     end
% end
% for i = 1:size(Lumped_NE_atOutcome_P, 1)
%     if ~isempty(find(Lumped_NE_atOutcome_P(i, :)<-0.6))
%         Lumped_NE_atOutcome_P(i, :) = NaN(1, size(Lumped_NE_atOutcome_P, 2));
%     end
% end
% for i = 1:size(Lumped_NE_atOutcome_R, 1)
%     if ~isempty(find(Lumped_NE_atOutcome_R(i, :)<-0.6))
%         Lumped_NE_atOutcome_R(i, :) = NaN(1, size(Lumped_NE_atOutcome_R, 2));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
subplot(2, 1, 1); hold on
ShadedPlot(TT2, nanmean(Lumped_NE_atOutcome_P, 1), [0 0 1], 1, SEM(Lumped_NE_atOutcome_P), [0.73 0.83 0.96])
plot(TT2, nanmean(Lumped_NE_atOutcome_P, 1), '--b', 'LineWidth', 1)
ShadedPlot(TT2, nanmean(Lumped_NE_atOutcome_R, 1), [0 0 1], 1, SEM(Lumped_NE_atOutcome_R), [0.73 0.83 0.96])
plot(TT2, nanmean(Lumped_NE_atOutcome_R, 1), '-b', 'LineWidth', 1)
ylabel('Z-score');
vline(0, '-k')
% figure; hold on
subplot(2, 1, 2); hold on
ShadedPlot(TT2, nanmean(Lumped_Ach_atOutcome_P, 1), [1 0 0], 1, SEM(Lumped_Ach_atOutcome_P), [0.9 0.8 0.7])
plot(TT2, nanmean(Lumped_Ach_atOutcome_P, 1), '--r', 'LineWidth', 1)
ShadedPlot(TT2, nanmean(Lumped_Ach_atOutcome_R, 1), [1 0 0], 1, SEM(Lumped_Ach_atOutcome_R), [0.9 0.8 0.7])
plot(TT2, nanmean(Lumped_Ach_atOutcome_R, 1), '-r', 'LineWidth', 1)
xlabel('Time after Outcome (s)')
ylabel('Z-score');
vline(0, '-k')
%%
% figure; hold on
subplot(1, 2, 1); 
ShadedPlot(TT2, mean([Lumped_NE_atOutcome_P; Lumped_NE_atOutcome_R], 1), [0 0 1], 1, SEM([Lumped_NE_atOutcome_P; Lumped_NE_atOutcome_R]), [0.73 0.83 0.96])
plot(TT2, mean([Lumped_NE_atOutcome_P; Lumped_NE_atOutcome_R], 1), '-b', 'LineWidth', 1)
ylabel('Z-score');
xlabel('Time after Outcome (s)')
vline(0, '-k')
subplot(1, 2, 2); 
ShadedPlot(TT2, mean([Lumped_Ach_atOutcome_P; Lumped_Ach_atOutcome_R], 1), [1 0 0], 1, SEM([Lumped_Ach_atOutcome_P; Lumped_Ach_atOutcome_R]), [0.9 0.8 0.7])
plot(TT2, mean([Lumped_Ach_atOutcome_P; Lumped_Ach_atOutcome_R], 1), '-r', 'LineWidth', 1)
xlabel('Time after Outcome (s)')
vline(0, '-k')
%% Plot the average amplitude of raw NE or ACh around the tone, before the outcomes
figure; hold on
% [WINst WINed] = [1. 2] for pre-onset mean; [2 4] for pre-outcome mean
WINst = 2*120+1;
WINed = 4*120;
MEAN_Ach_P = [];
MEAN_NE_P = [];
MEAN_Ach_R = [];
MEAN_NE_R = [];
% for subj = 1:length(subjIDX)
%     MEAN_Ach_P(end+1, 1) = nanmean(nanmean(Lumped_Ach_atOutcome_P(subjIDX{subj}, WINst:WINed), 2));
%     MEAN_NE_P(end+1, 1) = nanmean(nanmean(Lumped_NE_atOutcome_P(subjIDX{subj}, WINst:WINed), 2));
%     MEAN_Ach_R(end+1, 1) = nanmean(nanmean(Lumped_Ach_atOutcome_R(subjIDX{subj}, WINst:WINed), 2));
%     MEAN_NE_R(end+1, 1) = nanmean(nanmean(Lumped_NE_atOutcome_R(subjIDX{subj}, WINst:WINed), 2));
% end
% subplot(2, 1, 1); hold on
% bar(1:2, [nanmean(MEAN_Ach_P) nanmean(MEAN_Ach_R)])
% for i = 1:length(subjIDX)
%     plot(1:2, [MEAN_Ach_P(i) MEAN_Ach_R(i)], '-ok')
% end
% errorbar(1:2, [nanmean(MEAN_Ach_P) nanmean(MEAN_Ach_R)], [SEM(MEAN_Ach_P) SEM(MEAN_Ach_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Mean (z-score)')
% 
% subplot(2, 1, 2); hold on
% bar(1:2, [nanmean(MEAN_NE_P) nanmean(MEAN_NE_R)])
% for i = 1:length(subjIDX)
%     plot(1:2, [MEAN_NE_P(i) MEAN_NE_R(i)], '-ok')
% end
% errorbar(1:2, [nanmean(MEAN_NE_P) nanmean(MEAN_NE_R)], [SEM(MEAN_NE_P) SEM(MEAN_NE_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Mean (z-score)')

%%%%%%%%%% Session-based plot
subplot(2, 1, 1); hold on
bar(1:2, [nanmean(nanmean(Lumped_Ach_atOutcome_P(:, WINst:WINed), 2)) nanmean(nanmean(Lumped_Ach_atOutcome_R(:, WINst:WINed), 2))])
for i = 1:size(Lumped_Ach_atOutcome_P, 1)
    plot(1:2, [nanmean(Lumped_Ach_atOutcome_P(i, WINst:WINed), 2) nanmean(Lumped_Ach_atOutcome_R(i, WINst:WINed), 2)], '-ok')
end
errorbar(1:2, [nanmean(nanmean(Lumped_Ach_atOutcome_P(:, WINst:WINed), 2)) nanmean(nanmean(Lumped_Ach_atOutcome_R(:, WINst:WINed), 2))], ...
    [SEM(nanmean(Lumped_Ach_atOutcome_P(:, WINst:WINed), 2)) SEM(nanmean(Lumped_Ach_atOutcome_R(:, WINst:WINed), 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean (z-score)')

subplot(2, 1, 2); hold on
bar(1:2, [nanmean(nanmean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2)) nanmean(nanmean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2))])
for i = 1:size(Lumped_Ach_atOutcome_P, 1)
    plot(1:2, [nanmean(Lumped_NE_atOutcome_P(i, WINst:WINed), 2) nanmean(Lumped_NE_atOutcome_R(i, WINst:WINed), 2)], '-ok')
end
errorbar(1:2, [nanmean(nanmean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2)) nanmean(nanmean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2))], ...
    [SEM(nanmean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2)) SEM(nanmean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean (z-score)')
%% Find the trough time of raw traces punished trials (subject trough)
% looking for the trough in [-1s 0] -> outcome
WINst = 4*120+1;
WINed = 5*120;
TT = (-5*120+1:0*120)/120;
Trough_Ach_P = [];
Trough_NE_P = [];
tau_on_Ach = []
tau_off_Ach = [];
tau_on_NE = [];
tau_off_NE = [];
%%%%% Subject Troughs
for i = 1:length(subjIDX)
    [val, loc] = min(nanmean(Lumped_Ach_atOutcome_P(subjIDX{i}, WINst:WINed), 1));
    Trough_Ach_P(end+1, 1) = loc-1+WINst;
    [val, loc] = min(nanmean(Lumped_NE_atOutcome_P(subjIDX{i}, WINst:WINed), 1));
    Trough_NE_P(end+1, 1) = loc-1+WINst;
    %%% bio-exponential fit method
%     X = (1:length(TT))/120;
%     badfit_Ach = 1;
%     badfit_NE = 1;
%     while badfit_Ach
%         [x_fit, y_fit, tau_on, tau_off] = fit_bioexp(TT, mean(Lumped_Ach_atOutcome_P(subjIDX{i}, 1:length(TT)), 1));
%         [val, loc] = max(y_fit);
%         if X(loc)<=2 && X(loc)>0.01
%             badfit_Ach = 0;
%         end
%     end
%     Trough_Ach_P(end+1, 1) = -X(loc);
%     tau_on_Ach(end+1, 1) = tau_on;
%     tau_off_Ach(end+1, 1) = tau_off;
%     while badfit_NE
%         [x_fit, y_fit, tau_on, tau_off] = fit_bioexp(TT, mean(Lumped_NE_atOutcome_P(subjIDX{i}, 1:length(TT)), 1));
%         [val, loc] = max(y_fit);
%         if X(loc)<=2 && X(loc)>0.01
%             badfit_NE = 0;
%         end
%     end
%     Trough_NE_P(end+1, 1) = -X(loc);
%     tau_on_NE(end+1, 1) = tau_on;
%     tau_off_NE(end+1, 1) = tau_off;
end
%%%%% Session Troughs
% for i = 1:len
%     [val, loc] = min(nanmean(Lumped_Ach_atOutcome_P(i, WINst:WINed), 1));
%     Trough_Ach_P(end+1, 1) = loc-1+WINst;
%     [val, loc] = min(nanmean(Lumped_NE_atOutcome_P(i, WINst:WINed), 1));
%     Trough_NE_P(end+1, 1) = loc-1+WINst;
    %%%% bio-exponential fit method
%     X = (1:length(TT))/120;
%     badfit_Ach = 1;
%     badfit_NE = 1;
%     while badfit_Ach
%         [x_fit, y_fit, tau_on, tau_off] = fit_bioexp(TT, mean(Lumped_Ach_atOutcome_P(i, 1:length(TT)), 1));
%         [val, loc] = max(y_fit);
%         if X(loc)<3 && X(loc)>0.01
%             badfit_Ach = 0;
%         end
%     end
%     Trough_Ach_P(end+1, 1) = -X(loc);
%     tau_on_Ach(end+1, 1) = tau_on;
%     tau_off_Ach(end+1, 1) = tau_off;
%     while badfit_NE
%         [x_fit, y_fit, tau_on, tau_off] = fit_bioexp(TT, mean(Lumped_NE_atOutcome_P(i, 1:length(TT)), 1));
%         [val, loc] = max(y_fit);
%         if X(loc)<2 && X(loc)>0.01
%             badfit_NE = 0;
%         end
%     end
%     Trough_NE_P(end+1, 1) = -X(loc);    
%     tau_on_NE(end+1, 1) = tau_on;
%     tau_off_NE(end+1, 1) = tau_off;
% end
%%%%%%NOW GO BACK TO PREVIOUS SECTION AND PLOT %%%%%%%%%%%
% figure; hold on
% WINst = 1;
% MEAN_Ach_P = [];
% MEAN_NE_P = [];
% MEAN_Ach_R = [];
% MEAN_NE_R = [];
% for subj = 1:length(subjIDX)
%     WINed = Trough_Ach_P(subj);
%     MEAN_Ach_P(end+1, 1) = nanmean(nanmean(Lumped_Ach_atOutcome_P(subjIDX{subj}, WINst:WINed), 2));
%     MEAN_Ach_R(end+1, 1) = nanmean(nanmean(Lumped_Ach_atOutcome_R(subjIDX{subj}, WINst:WINed), 2));
%     WINed = Trough_NE_P(subj);
%     MEAN_NE_P(end+1, 1) = nanmean(nanmean(Lumped_NE_atOutcome_P(subjIDX{subj}, WINst:WINed), 2));
%     MEAN_NE_R(end+1, 1) = nanmean(nanmean(Lumped_NE_atOutcome_R(subjIDX{subj}, WINst:WINed), 2));
% end
% subplot(2, 1, 1); hold on
% bar(1:2, [nanmean(MEAN_Ach_P) nanmean(MEAN_Ach_R)])
% for i = 1:length(subjIDX)
%     plot(1:2, [MEAN_Ach_P(i) MEAN_Ach_R(i)], '-ok')
% end
% errorbar(1:2, [nanmean(MEAN_Ach_P) nanmean(MEAN_Ach_R)], [SEM(MEAN_Ach_P) SEM(MEAN_Ach_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Mean (z-score)')
% 
% subplot(2, 1, 2); hold on
% bar(1:2, [nanmean(MEAN_NE_P) nanmean(MEAN_NE_R)])
% for i = 1:length(subjIDX)
%     plot(1:2, [MEAN_NE_P(i) MEAN_NE_R(i)], '-ok')
% end
% errorbar(1:2, [nanmean(MEAN_NE_P) nanmean(MEAN_NE_R)], [SEM(MEAN_NE_P) SEM(MEAN_NE_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Mean (z-score)')
%% Store Trough for Saline
Trough_NE_Saline = Trough_NE_P;
Trough_Ach_Saline = Trough_Ach_P;
%% Store Trough for CNO
Trough_NE_CNO = Trough_NE_P;
Trough_Ach_CNO = Trough_Ach_P;
%% Store tau of bi-exponential fitting for Saline
tau_on_NE_Saline = tau_on_NE;
tau_off_NE_Saline = tau_off_NE;
tau_on_Ach_Saline = tau_on_Ach;
tau_off_Ach_Saline = tau_off_Ach;
%% Store tau of bi-exponential fitting for CNO
tau_on_NE_CNO = tau_on_NE;
tau_off_NE_CNO = tau_off_NE;
tau_on_Ach_CNO = tau_on_Ach;
tau_off_Ach_CNO = tau_off_Ach;
%% Plot bi exponential tau of failed for Saline vs CNO
%%% Bi-exponential fitting method plot
figure; hold on
bar(1:2, [mean(tau_on_NE_Saline) mean(tau_on_NE_CNO)])
bar(4:5, [mean(tau_on_Ach_Saline) mean(tau_on_Ach_CNO)])
% for i = 1:length(subjIDX)
%     plot(1:2, [Trough_NE_Saline(i) Trough_NE_CNO(i)], '-ok')
%     plot(4:5, [Trough_Ach_Saline(i) Trough_Ach_CNO(i)], '-ok')
% end
scatter(ones(length(tau_on_NE_Saline), 1), tau_on_NE_Saline)
scatter(2*ones(length(tau_on_NE_CNO), 1), tau_on_NE_CNO)
scatter(4*ones(length(tau_on_Ach_Saline), 1), tau_on_Ach_Saline)
scatter(5*ones(length(tau_on_Ach_CNO), 1), tau_on_Ach_CNO)
errorbar([1:2 4:5], [mean(tau_on_NE_Saline) mean(tau_on_NE_CNO) mean(tau_on_Ach_Saline) mean(tau_on_Ach_CNO)], ...
    [SEM(tau_on_NE_Saline) SEM(tau_on_NE_CNO) SEM(tau_on_Ach_Saline) SEM(tau_on_Ach_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
bar(1:2, -[mean(tau_off_NE_Saline) mean(tau_off_NE_CNO)])
bar(4:5, -[mean(tau_off_Ach_Saline) mean(tau_off_Ach_CNO)])
% for i = 1:length(subjIDX)
%     plot(1:2, [Trough_NE_Saline(i) Trough_NE_CNO(i)], '-ok')
%     plot(4:5, [Trough_Ach_Saline(i) Trough_Ach_CNO(i)], '-ok')
% end
scatter(ones(length(tau_off_NE_Saline), 1), -tau_off_NE_Saline)
scatter(2*ones(length(tau_off_NE_CNO), 1), -tau_off_NE_CNO)
scatter(4*ones(length(tau_off_Ach_Saline), 1), -tau_off_Ach_Saline)
scatter(5*ones(length(tau_off_Ach_CNO), 1), -tau_off_Ach_CNO)
errorbar([1:2 4:5], -[mean(tau_off_NE_Saline) mean(tau_off_NE_CNO) mean(tau_off_Ach_Saline) mean(tau_off_Ach_CNO)], ...
    [SEM(tau_off_NE_Saline) SEM(tau_off_NE_CNO) SEM(tau_off_Ach_Saline) SEM(tau_off_Ach_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
%% Plot trough time of failed for Saline vs CNO
%%% Bi-exponential fitting method plot
% figure; hold on
% bar(1:2, [mean(Trough_NE_Saline) mean(Trough_NE_CNO)])
% bar(4:5, [mean(Trough_Ach_Saline) mean(Trough_Ach_CNO)])
% % for i = 1:length(subjIDX)
% %     plot(1:2, [Trough_NE_Saline(i) Trough_NE_CNO(i)], '-ok')
% %     plot(4:5, [Trough_Ach_Saline(i) Trough_Ach_CNO(i)], '-ok')
% % end
% scatter(ones(length(Trough_NE_Saline), 1), Trough_NE_Saline)
% scatter(2*ones(length(Trough_NE_CNO), 1), Trough_NE_CNO)
% scatter(4*ones(length(Trough_Ach_Saline), 1), Trough_Ach_Saline)
% scatter(5*ones(length(Trough_Ach_CNO), 1), Trough_Ach_CNO)
% errorbar([1:2 4:5], [mean(Trough_NE_Saline) mean(Trough_NE_CNO) mean(Trough_Ach_Saline) mean(Trough_Ach_CNO)], ...
%     [SEM(Trough_NE_Saline) SEM(Trough_NE_CNO) SEM(Trough_Ach_Saline) SEM(Trough_Ach_CNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

%%%% no fitting
figure; hold on
bar(1, -(5*120-mean(Trough_NE_Saline))/120)
bar(2, -(5*120-mean(Trough_NE_CNO))/120)
bar(4, -(5*120-mean(Trough_Ach_Saline))/120)
bar(5, -(5*120-mean(Trough_Ach_CNO))/120)
% for i = 1:length(subjIDX)
%     plot(1:2, -(5*120-[Trough_NE_Saline(i) Trough_NE_CNO(i)])/120, '-ok')
%     plot(4:5, -(5*120-[Trough_Ach_Saline(i) Trough_Ach_CNO(i)])/120, '-ok')
% end
scatter(ones(length(Trough_NE_Saline), 1), -(5*120-Trough_NE_Saline)/120)
scatter(2*ones(length(Trough_NE_CNO), 1), -(5*120-Trough_NE_CNO)/120)
scatter(4*ones(length(Trough_Ach_Saline), 1), -(5*120-Trough_Ach_Saline)/120)
scatter(5*ones(length(Trough_Ach_CNO), 1), -(5*120-Trough_Ach_CNO)/120)
errorbar([1:2 4:5], -(5*120-[mean(Trough_NE_Saline) mean(Trough_NE_CNO) mean(Trough_Ach_Saline) mean(Trough_Ach_CNO)])/120, ...
    [SEM(Trough_NE_Saline) SEM(Trough_NE_CNO) SEM(Trough_Ach_Saline) SEM(Trough_Ach_CNO)]/120, 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Trough time (s)')
%% Plot the trough time of failed for all WT and Saline
figure; hold on
bar(1, -(5*120-mean(Trough_NE_P))/120)
bar(2, -(5*120-mean(Trough_Ach_P))/120)
for i = 1:length(subjIDX)
    plot(1:2, -(5*120-[Trough_NE_P(i) Trough_Ach_P(i)])/120, '-ok')
end
errorbar(1:2, -(5*120-[mean(Trough_NE_P) mean(Trough_Ach_P)])/120, [SEM(Trough_NE_P) SEM(Trough_Ach_P)]/120, 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Trough time (s)')
%% Find the trough time of raw traces successful trials (subject trough)
% looking for the trough in [0 1s] -> outcome
WINst = 5*120;
WINed = 6*120+1;
Trough_Ach_R = [];
Trough_NE_R = [];
for i = 1:length(subjIDX)
    [val, loc] = min(nanmean(Lumped_Ach_atOutcome_R(subjIDX{i}, WINst:WINed), 1));
    Trough_Ach_R(end+1, 1) = loc-1+WINst;
    [val, loc] = min(nanmean(Lumped_NE_atOutcome_R(subjIDX{i}, WINst:WINed), 1));
    Trough_NE_R(end+1, 1) = loc-1+WINst;
end
%% Plot the trough time of successful
figure; hold on
bar(1, -(5*120-mean(Trough_NE_R))/120)
bar(2, -(5*120-mean(Trough_Ach_R))/120)
for i = 1:length(subjIDX)
    plot(1:2, -(5*120-[Trough_NE_R(i) Trough_Ach_R(i)])/120, '-ok')
end
errorbar(1:2, -(5*120-[mean(Trough_NE_R) mean(Trough_Ach_R)])/120, [SEM(Trough_NE_R) SEM(Trough_Ach_R)]/120, 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Trough time (s)')
%% Find the slope of raw sensor traces at outcomes
[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
    MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 0, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 0, 1);

Lumped_NE_atOutcome_P = NaN(len, 841); % Punished but withheld >4s NE
Lumped_Ach_atOutcome_P = NaN(len, 841); % Punished but withheld >4s Ach
Lumped_NE_atOutcome_R = NaN(len, 841); % Rewarded NE
Lumped_Ach_atOutcome_R = NaN(len, 841); % Rewarded Ach

tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 4;
bandSensor = [0.3 1.2];
WIN_spikecnt = 4;

cnt = 0;
for i = OFCIDX(:).'
    cnt = cnt+1;
    z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, [0.1 10], 'bandpass'));
%     z_NE = runningmean(e1.MetaData.NE_470{i}, 0, 120, 30, 1);
    z_Ach = zscore(Filter(e1.MetaData.Ach_470{i}, 120, 2, [0.1 10], 'bandpass'));
%     z_Ach = runningmean(e1.MetaData.Ach_470{i}, 0, 120, 30, 1);
%     z_Ach = e1.MetaData.Ach_F{i};

    i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass'); i_NE = zscore(i_NE(e1.MetaData.Pupil_in_sensor{i}));
    i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass'); i_Ach = zscore(i_Ach(e1.MetaData.Pupil_in_sensor{i}));

    m = load(MetaData_files{i});
    mx = load(MetaDataX_files{i});

    session_NE_atOutcome_P = NaN(height(m.MetaData), 841);
    session_Ach_atOutcome_P = NaN(height(m.MetaData), 841);
    session_NE_atOutcome_R = NaN(height(m.MetaData), 841);
    session_Ach_atOutcome_R = NaN(height(m.MetaData), 841);
    for k = 1:height(m.MetaData)
        stTone_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, m.MetaData.Tone_Onset(k));
        if m.MetaData.Punish_Onset(k)-m.MetaData.Tone_Onset(k)>4
            if m.MetaData.Punish_Onset(k)+5<=e1.MetaData.Sensor_time{i}(end)
                stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, m.MetaData.Punish_Onset(k));
                session_NE_atOutcome_P(k, :) = z_NE(stOutcome_sensor-720:stOutcome_sensor+120);
                session_Ach_atOutcome_P(k, :) = z_Ach(stOutcome_sensor-720:stOutcome_sensor+120);
            end
        end
        if ~isnan(m.MetaData.Reward_Onset(k)) 
            if m.MetaData.Reward_Onset(k)+5<=e1.MetaData.Sensor_time{i}(end)
                stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, m.MetaData.Reward_Onset(k));
                session_NE_atOutcome_R(k, :) = z_NE(stOutcome_sensor-720:stOutcome_sensor+120);
                session_Ach_atOutcome_R(k, :) = z_Ach(stOutcome_sensor-720:stOutcome_sensor+120);
            end
        end
    end
    Lumped_NE_atOutcome_P(cnt, :) = nanmean(session_NE_atOutcome_P, 1); 
    Lumped_Ach_atOutcome_P(cnt, :) = nanmean(session_Ach_atOutcome_P, 1); 
    Lumped_NE_atOutcome_R(cnt, :) = nanmean(session_NE_atOutcome_R, 1); 
    Lumped_Ach_atOutcome_R(cnt, :) = nanmean(session_Ach_atOutcome_R, 1);
end
toc
%% Find the slope of raw traces at outcomes
% connecting the slope between [-5s -4s] and [-2s -1s] -> outcome; generating the [-5 -1] slope
% connecting the slope between [-3.5s -2.5s] and [-1.5s -0.5s] -> outcome; generating the [-3 -1] slope
WINst1 = 2.5*120;
WINed1 = 3.5*120;
WINst2 = 4.5*120;
WINed2 = 5.5*120;
Slope_Ach_P = [];
Slope_NE_P = [];
Slope_Ach_R = [];
Slope_NE_R = [];
for i = 1:length(subjIDX)
    trace_Ach_P = medfilt1(nanmean(Lumped_Ach_atOutcome_P(subjIDX{i}, :), 1), 40);
    Slope_Ach_P(end+1, 1) = (mean(trace_Ach_P(WINst1:WINed1))-mean(trace_Ach_P(WINst2:WINed2)))/((WINed1-WINed2)/120);
    trace_Ach_R = medfilt1(nanmean(Lumped_Ach_atOutcome_R(subjIDX{i}, :), 1), 40);
    Slope_Ach_R(end+1, 1) = (mean(trace_Ach_R(WINst1:WINed1))-mean(trace_Ach_R(WINst2:WINed2)))/((WINed1-WINed2)/120);
    trace_NE_P = medfilt1(nanmean(Lumped_NE_atOutcome_P(subjIDX{i}, :), 1), 40);
    Slope_NE_P(end+1, 1) = (mean(trace_NE_P(WINst1:WINed1))-mean(trace_NE_P(WINst2:WINed2)))/((WINed1-WINed2)/120);
    trace_NE_R = medfilt1(nanmean(Lumped_NE_atOutcome_R(subjIDX{i}, :), 1), 40);
    Slope_NE_R(end+1, 1) = (mean(trace_NE_R(WINst1:WINed1))-mean(trace_NE_R(WINst2:WINed2)))/((WINed1-WINed2)/120);
end

%%%%%%%%%%% Session-based plot
% sessSlope_Ach_P = [];
% sessSlope_NE_P = [];
% sessSlope_Ach_R = [];
% sessSlope_NE_R = [];
% for i = 1:len
%     trace_Ach_P = medfilt1(Lumped_Ach_atOutcome_P(i, :), 40);
%     sessSlope_Ach_P(end+1, 1) = (mean(trace_Ach_P(WINst1:WINed1))-mean(trace_Ach_P(WINst2:WINed2)))/((WINed1-WINed2)/120);
%     trace_Ach_R = medfilt1(Lumped_Ach_atOutcome_R(i, :), 40);
%     sessSlope_Ach_R(end+1, 1) = (mean(trace_Ach_R(WINst1:WINed1))-mean(trace_Ach_R(WINst2:WINed2)))/((WINed1-WINed2)/120);
%     trace_NE_P = medfilt1(Lumped_NE_atOutcome_P(i, :), 40);
%     sessSlope_NE_P(end+1, 1) = (mean(trace_NE_P(WINst1:WINed1))-mean(trace_NE_P(WINst2:WINed2)))/((WINed1-WINed2)/120);
%     trace_NE_R = medfilt1(Lumped_NE_atOutcome_R(i, :), 40);
%     sessSlope_NE_R(end+1, 1) = (mean(trace_NE_R(WINst1:WINed1))-mean(trace_NE_R(WINst2:WINed2)))/((WINed1-WINed2)/120);
% end
%% Plot the slope 
figure; hold on
subplot(2, 1, 1); hold on
bar(1:2, [mean(Slope_Ach_P) mean(Slope_Ach_R)])
for i = 1:length(subjIDX)
    plot(1:2, [Slope_Ach_P(i) Slope_Ach_R(i)], '-ok')
end
errorbar(1:2, [mean(Slope_Ach_P) mean(Slope_Ach_R)], [SEM(Slope_Ach_P) SEM(Slope_Ach_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Slope (a.u.)')

subplot(2, 1, 2); hold on
bar(1:2, [mean(Slope_NE_P) mean(Slope_NE_R)])
for i = 1:length(subjIDX)
    plot(1:2, [Slope_NE_P(i) Slope_NE_R(i)], '-ok')
end
errorbar(1:2, [mean(Slope_NE_P) mean(Slope_NE_R)], [SEM(Slope_NE_P) SEM(Slope_NE_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Slope (a.u.)')

%%%%%%%%%%% Session-based plot
% subplot(2, 1, 1); hold on
% bar(1:2, [mean(sessSlope_Ach_P) mean(sessSlope_Ach_R)])
% for i = 1:size(sessSlope_Ach_P, 1)
%     plot(1:2, [sessSlope_Ach_P(i) sessSlope_Ach_R(i)], '-ok')
% end
% errorbar(1:2, [mean(sessSlope_Ach_P) mean(sessSlope_Ach_R)], [SEM(sessSlope_Ach_P) SEM(sessSlope_Ach_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Slope (a.u.)')
% 
% subplot(2, 1, 2); hold on
% bar(1:2, [mean(sessSlope_NE_P) mean(sessSlope_NE_R)])
% for i = 1:size(Lumped_Ach_atOutcome_P, 1)
%     plot(1:2, [sessSlope_NE_P(i) sessSlope_NE_R(i)], '-ok')
% end
% errorbar(1:2, [mean(sessSlope_NE_P) mean(sessSlope_NE_R)], [SEM(sessSlope_NE_P) SEM(sessSlope_NE_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
% ylabel('Slope (a.u.)')
%% Find the peak value and peak time of raw traces failed trials around tone (subject peak)
% looking for the peak in [0 1s] -> tone onset
% MAKE SURE THE MEAN_Ach_P IS THE PRE-TONE MEAN FIRST!!!!
WINst = 2*120;
WINed = 3*120+1;
Peak_Ach_P = [];
Peak_NE_P = [];
Peak_Ach_R = [];
Peak_NE_R = [];
for i = 1:length(subjIDX)
    [val, loc] = max(medfilt1(nanmean(Lumped_Ach_atTone_P(subjIDX{i}, WINst:WINed), 1), 4));
%     [val, loc] = max(nanmean(Lumped_Ach_atTone_P(subjIDX{i}, WINst:WINed), 1));
    Peak_Ach_P(end+1, :) = [loc/120 val-MEAN_Ach_P(i)];
%     [val, loc] = max(nanmean(Lumped_Ach_atTone_R(subjIDX{i}, WINst:WINed), 1));
    [val, loc] = max(medfilt1(nanmean(Lumped_Ach_atTone_R(subjIDX{i}, WINst:WINed), 1), 4));
    Peak_Ach_R(end+1, :) = [loc/120 val-MEAN_Ach_R(i)];
%     [val, loc] = max(nanmean(Lumped_NE_atTone_P(subjIDX{i}, WINst:WINed), 1));
    [val, loc] = max(medfilt1(nanmean(Lumped_NE_atTone_P(subjIDX{i}, WINst:WINed), 1), 4));
    Peak_NE_P(end+1, :) = [loc/120 val-MEAN_NE_P(i)];
%     [val, loc] = max(nanmean(Lumped_NE_atTone_R(subjIDX{i}, WINst:WINed), 1));
    [val, loc] = max(medfilt1(nanmean(Lumped_NE_atTone_R(subjIDX{i}, WINst:WINed), 1), 4));
    Peak_NE_R(end+1, :) = [loc/120 val-MEAN_NE_R(i)];
end
%% Plot the peak value
figure; hold on
subplot(2, 1, 1); hold on
bar(1:2, [mean(Peak_Ach_P(:, 2)) mean(Peak_Ach_R(:, 2))])
for i = 1:length(subjIDX)
    plot(1:2, [Peak_Ach_P(i, 2) Peak_Ach_R(i, 2)], '-ok')
end
errorbar(1:2, [mean(Peak_Ach_P(:, 2)) mean(Peak_Ach_R(:, 2))], [SEM(Peak_Ach_P(:, 2)) SEM(Peak_Ach_R(:, 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Post-onset peak (z-score)')

subplot(2, 1, 2); hold on
bar(1:2, [mean(Peak_NE_P(:, 2)) mean(Peak_NE_R(:, 2))])
for i = 1:length(subjIDX)
    plot(1:2, [Peak_NE_P(i, 2) Peak_NE_R(i, 2)], '-ok')
end
errorbar(1:2, [mean(Peak_NE_P(:, 2)) mean(Peak_NE_R(:, 2))], [SEM(Peak_NE_P(:, 2)) SEM(Peak_NE_R(:, 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Post-onset peak (z-score)')
%% Plot the peak latency
figure; hold on
subplot(2, 1, 1); hold on
bar(1:2, [mean(Peak_Ach_P(:, 1)) mean(Peak_Ach_R(:, 1))])
for i = 1:length(subjIDX)
    plot(1:2, [Peak_Ach_P(i, 1) Peak_Ach_R(i, 1)], '-ok')
end
errorbar(1:2, [mean(Peak_Ach_P(:, 1)) mean(Peak_Ach_R(:, 1))], [SEM(Peak_Ach_P(:, 1)) SEM(Peak_Ach_R(:, 1))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Post-onset peak (z-score)')

subplot(2, 1, 2); hold on
bar(1:2, [mean(Peak_NE_P(:, 1)) mean(Peak_NE_R(:, 1))])
for i = 1:length(subjIDX)
    plot(1:2, [Peak_NE_P(i, 1) Peak_NE_R(i, 1)], '-ok')
end
errorbar(1:2, [mean(Peak_NE_P(:, 1)) mean(Peak_NE_R(:, 1))], [SEM(Peak_NE_P(:, 1)) SEM(Peak_NE_R(:, 1))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Post-onset peak (z-score)')
%% Lumped raw sensor traces before outcome under Saline
MEAN_NE_noCNO = [];
MEAN_Ach_noCNO = [];
for i = 1:length(subjIDX)
    MEAN_NE_noCNO(end+1, :) = nanmean([subj_NE_atOutcome_P(i, :); subj_NE_atOutcome_R(i, :)], 1);
    MEAN_Ach_noCNO(end+1, :) = nanmean([subj_Ach_atOutcome_P(i, :); subj_Ach_atOutcome_R(i, :)], 1);
end


%% Lumped raw sensor traces before outcome under CNO
MEAN_NE_CNO = [];
MEAN_Ach_CNO = [];
for i = 1:length(subjIDX)
    MEAN_NE_CNO(end+1, :) = nanmean([subj_NE_atOutcome_P(i, :); subj_NE_atOutcome_R(i, :)], 1);
    MEAN_Ach_CNO(end+1, :) = nanmean([subj_Ach_atOutcome_P(i, :); subj_Ach_atOutcome_R(i, :)], 1);
end

%% Plot raw sensor traces before outcome Saline vs CNO
figure; hold on
subplot(2, 1, 1); hold on
ShadedPlot(TT2, mean(MEAN_NE_noCNO, 1), [0 0 1], 1, SEM(MEAN_NE_noCNO), [0.73 0.83 0.96])
ShadedPlot(TT2, mean(MEAN_NE_CNO, 1), [0 0 1], 1, SEM(MEAN_NE_CNO), [0.73 0.83 0.96])
plot(TT2, mean(MEAN_NE_noCNO, 1), '-b', 'LineWidth', 1)
plot(TT2, mean(MEAN_NE_CNO, 1), '-c', 'LineWidth', 1)
ylabel('Z-score');
vline(0, '-k')

subplot(2, 1, 2); hold on
ShadedPlot(TT2, mean(MEAN_Ach_noCNO, 1), [0 0 1], 1, SEM(MEAN_Ach_noCNO), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(MEAN_Ach_CNO, 1), [0 0 1], 1, SEM(MEAN_Ach_CNO), [0.9 0.8 0.7])
plot(TT2, mean(MEAN_Ach_noCNO, 1), '-r', 'LineWidth', 1)
plot(TT2, mean(MEAN_Ach_CNO, 1), '-m', 'LineWidth', 1)
ylabel('Z-score');
xlabel('Time after Outcome (s)')
vline(0, '-k')
%% Mean value of raw trace ratio between failed and successful trials
WINst = 3.5*120+1;
WINed = 4*120;
% % MEAN_NE_ratio = [];
% % MEAN_Ach_ratio = [];
% MEAN_NE_diff = [];
% MEAN_Ach_diff = [];
% for i = 1:length(subjIDX)
% %     MEAN_NE_ratio(end+1, 1) = mean(mean(Lumped_NE_atOutcome_R(subjIDX{i}, WINst:WINed), 2)./mean(Lumped_NE_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
%     MEAN_NE_diff(end+1, 1) = mean(mean(Lumped_NE_atOutcome_R(subjIDX{i}, WINst:WINed), 2)-mean(Lumped_NE_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
% %     MEAN_Ach_ratio(end+1, 1) = mean(mean(Lumped_Ach_atOutcome_R(subjIDX{i}, WINst:WINed), 2)./mean(Lumped_Ach_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
%     MEAN_Ach_diff(end+1, 1) = mean(mean(Lumped_Ach_atOutcome_R(subjIDX{i}, WINst:WINed), 2)-mean(Lumped_Ach_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
% end

MEAN_NE_ratio = mean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2)./mean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2);
MEAN_Ach_ratio = mean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2)./mean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2);
MEAN_NE_diff = mean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2)-mean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2);
MEAN_Ach_diff = mean(Lumped_Ach_atOutcome_R(:, WINst:WINed), 2)-mean(Lumped_Ach_atOutcome_P(:, WINst:WINed), 2);

%% Mean value of raw trace ratio between failed and successful trials
WINst = 3.5*120+1;
WINed = 4*120;
% MEAN_NE_ratioCNO = [];
% MEAN_Ach_ratioCNO = [];
% % MEAN_NE_diffCNO = [];
% % MEAN_Ach_diffCNO = [];
% for i = 1:length(subjIDX)
%     MEAN_NE_ratioCNO(end+1, 1) = mean(mean(Lumped_NE_atOutcome_R(subjIDX{i}, WINst:WINed), 2)./mean(Lumped_NE_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
% %     MEAN_NE_diffCNO(end+1, 1) = mean(mean(Lumped_NE_atOutcome_R(subjIDX{i}, WINst:WINed), 2)-mean(Lumped_NE_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
%     MEAN_Ach_ratioCNO(end+1, 1) = mean(mean(Lumped_Ach_atOutcome_R(subjIDX{i}, WINst:WINed), 2)./mean(Lumped_Ach_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
% %     MEAN_Ach_diffCNO(end+1, 1) = mean(mean(Lumped_Ach_atOutcome_R(subjIDX{i}, WINst:WINed), 2)-mean(Lumped_Ach_atOutcome_P(subjIDX{i}, WINst:WINed), 2));
% end

MEAN_NE_ratioCNO = mean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2)./mean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2);
MEAN_Ach_ratioCNO = mean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2)./mean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2);
MEAN_NE_diffCNO = mean(Lumped_NE_atOutcome_R(:, WINst:WINed), 2)-mean(Lumped_NE_atOutcome_P(:, WINst:WINed), 2);
MEAN_Ach_diffCNO = mean(Lumped_Ach_atOutcome_R(:, WINst:WINed), 2)-mean(Lumped_Ach_atOutcome_P(:, WINst:WINed), 2);
%% Plot raw trace mean ratio
figure; hold on
bar(1, mean(MEAN_NE_ratio))
bar(2, mean(MEAN_NE_ratioCNO));
bar(4, mean(MEAN_Ach_ratio));
bar(5, mean(MEAN_Ach_ratioCNO));
for i = 1:length(MEAN_NE_ratio)
%     plot(1:2, [MEAN_NE_ratio(i) MEAN_NE_ratioCNO(i)], '-ok')
%     plot(4:5, [MEAN_Ach_ratio(i) MEAN_Ach_ratioCNO(i)], '-ok')
end
scatter(ones(length(MEAN_NE_ratio), 1), MEAN_NE_ratio, 8, 'MarkerEdgeColor', 'k')
scatter(2*ones(length(MEAN_NE_ratioCNO), 1), MEAN_NE_ratioCNO, 8, 'MarkerEdgeColor', 'k')
scatter(4*ones(length(MEAN_Ach_ratio), 1), MEAN_Ach_ratio, 8, 'MarkerEdgeColor', 'k')
scatter(5*ones(length(MEAN_Ach_ratioCNO), 1), MEAN_Ach_ratioCNO, 8, 'MarkerEdgeColor', 'k')
errorbar([1 2 4 5], [mean(MEAN_NE_ratio) mean(MEAN_NE_ratioCNO) mean(MEAN_Ach_ratio) mean(MEAN_Ach_ratioCNO)], ...
    [SEM(MEAN_NE_ratio) SEM(MEAN_NE_ratioCNO) SEM(MEAN_Ach_ratio) SEM(MEAN_Ach_ratioCNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

figure; hold on
bar(1, mean(MEAN_NE_diff))
bar(2, mean(MEAN_NE_diffCNO));
bar(4, mean(MEAN_Ach_diff));
bar(5, mean(MEAN_Ach_diffCNO));
for i = 1:length(MEAN_NE_diff)
%     plot(1:2, [MEAN_NE_diff(i) MEAN_NE_diffCNO(i)], '-ok')
%     plot(4:5, [MEAN_Ach_diff(i) MEAN_Ach_diffCNO(i)], '-ok')
end
scatter(ones(length(MEAN_NE_diff), 1), MEAN_NE_diff, 8, 'MarkerEdgeColor', 'k')
scatter(2*ones(length(MEAN_NE_diffCNO), 1), MEAN_NE_diffCNO, 8, 'MarkerEdgeColor', 'k')
scatter(4*ones(length(MEAN_Ach_diff), 1), MEAN_Ach_diff, 8, 'MarkerEdgeColor', 'k')
scatter(5*ones(length(MEAN_Ach_diffCNO), 1), MEAN_Ach_diffCNO, 8, 'MarkerEdgeColor', 'k')
errorbar([1 2 4 5], [mean(MEAN_NE_diff) mean(MEAN_NE_diffCNO) mean(MEAN_Ach_diff) mean(MEAN_Ach_diffCNO)], ...
    [SEM(MEAN_NE_diff) SEM(MEAN_NE_diffCNO) SEM(MEAN_Ach_diff) SEM(MEAN_Ach_diffCNO)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);









