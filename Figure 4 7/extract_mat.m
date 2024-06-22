[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);
session_NE = cell(len, 1);
session_Ach = cell(len, 1);

session_sync = cell(len, 1);
% session_sync_RC1 = cell(len, 1);
% session_sync_RC2 = cell(len, 1);
% session_sync_RC3 = cell(len, 1);

session_spiketime = cell(len, 1);
session_spiketime_P = cell(len, 1);
session_spiketime_R = cell(len, 1);
session_OutcomeIDX_P = cell(len, 1);
session_OutcomeIDX_R = cell(len, 1);
session_sensorIDX_shuffle = cell(len, 1);
session_spike2P = cell(len, 1); % stores the spike times relative to the outcome onset
session_spike2R = cell(len, 1); % stores the spike times relative to the outcome onset
session_Nspike_P = cell(len, 1); % stores the number of spikes in the period before punishment
session_Nspike_R = cell(len, 1); % stores the number of spikes in the period before success
session_SHFspike_P = cell(len, 1);
session_SHFspike_R = cell(len, 1);
session_spike_TonePeriod_P = cell(len, 1);
session_spike_TonePeriod_R = cell(len, 1);
session_spike_Trial2Outcome_P = cell(len, 1);
session_spike_Trial2Outcome_R = cell(len, 1);

session_NE_atOutcome_P = cell(len, 1); % Punished but withheld >4s NE
session_Ach_atOutcome_P = cell(len, 1); % Punished but withheld >4s Ach
session_NA_atOutcome_P = cell(len, 1);
session_DiffNA_atOutcome_P = cell(len, 1);
session_NE_atOutcome_R = cell(len, 1); % Rewarded NE
session_Ach_atOutcome_R = cell(len, 1); % Rewarded Ach
session_NA_atOutcome_R = cell(len, 1);
session_DiffNA_atOutcome_R = cell(len, 1);

session_sync_atOutcome_NS = cell(len, 1);
session_sync_Free_NS = cell(len, 1);

session_sync_LDA_P = cell(len, 1); % LDA transformed features
session_sync_LDA_R = cell(len, 1); % LDA transformed features

session_NEAch_spkcnt_atOutcome_P = cell(len, 1);
session_NEAch_spkcnt_atOutcome_R = cell(len, 1);
session_NEAch_sync_atOutcome_P = cell(len, 1);
session_NEAch_sync_atOutcome_R = cell(len, 1);
session_NEAch_sync_TonePeriod_P = cell(len, 1);
session_NEAch_sync_TonePeriod_R = cell(len, 1);
session_NEAch_sync_Trial2Outcome_P = cell(len, 1);
session_NEAch_sync_Trial2Outcome_R = cell(len, 1);
session_NEAch_FR_atOutcome_P = cell(len, 1);
session_NEAch_FR_atOutcome_R = cell(len, 1);

session_NEAch_syncRC1_TonePeriod_P = cell(len, 1);
session_NEAch_syncRC1_TonePeriod_R = cell(len, 1);
session_NEAch_syncRC2_TonePeriod_P = cell(len, 1);
session_NEAch_syncRC2_TonePeriod_R = cell(len, 1);
session_NEAch_syncRC3_TonePeriod_P = cell(len, 1);
session_NEAch_syncRC3_TonePeriod_R = cell(len, 1);

session_NEAch_syncENV_TonePeriod_P = cell(len, 1);
session_NEAch_syncENV_TonePeriod_R = cell(len, 1);

session_NE_TonePeriod_P = cell(len, 1);
session_NE_TonePeriod_R = cell(len, 1);
session_Ach_TonePeriod_P = cell(len, 1);
session_Ach_TonePeriod_R = cell(len, 1);

session_NE_LickFreePeriod_P = cell(len, 1);
session_NE_LickFreePeriod_R = cell(len, 1);
session_Ach_LickFreePeriod_P =  cell(len, 1);
session_Ach_LickFreePeriod_R = cell(len, 1);
session_NEAch_sync_LickFreePeriod_P = cell(len, 1);
session_NEAch_sync_LickFreePeriod_R = cell(len, 1);

winsz = 5;
tmin = -200; tmax = 1000;
tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 3; 
bandSensor = [0.4 0.8];
WINspkcnt_st = -5;
WINspkcnt_ed = 0;
WINfree_st = -5;
WINfree_ed = 0;
EDGE_FR = -3.5:0.5:-1;
num_SHF = 1;
fs = 120;
cnt = 0;

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
%         z_NE = zscore(e1.MetaData.NE_470{i});
%         z_Ach = zscore(e1.MetaData.Ach_470{i});
%     i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass'); i_NE = Filter(i_NE(e1.MetaData.Pupil_in_sensor{i}), 10, 4, 3.5, 'low');
    
%     i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, 0.3, 'high'); i_NE = Filter(i_NE, 120, 2, 1.2, 'low'); i_NE = i_NE(S2P_IDX);
%     i_NE = Filter(e1.MetaData.NE_470{i}, 120, 3, bandSensor, 'bandpass');
%     i_NE = Filter(e1.MetaData.NE_F{i}, 120, 2, bandSensor, 'bandpass'); i_NE = Filter(i_NE(e1.MetaData.Pupil_in_sensor{i}), 10, 4, 3.5, 'low');
%     i_NE = Filter(e1.MetaData.NE_470{i}(S2P_IDX), 10, 4, 3.5, 'low'); i_NE = zscore(Filter(i_NE, 10, 4, bandSensor, 'bandpass'));
%     i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass'); i_Ach = Filter(i_Ach(e1.MetaData.Pupil_in_sensor{i}), 10, 4, 3.5, 'low');
    
%     i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, 0.3, 'high'); i_Ach = Filter(i_Ach, 120, 2, 1.2, 'low'); i_Ach = i_Ach(S2P_IDX);
%     i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 3, bandSensor, 'bandpass');
%     i_Ach = Filter(e1.MetaData.Ach_F{i}, 120, 2, bandSensor, 'bandpass'); i_Ach = Filter(i_Ach(e1.MetaData.Pupil_in_sensor{i}), 10, 4, 3.5, 'low');
%     i_Ach = Filter(e1.MetaData.Ach_470{i}(S2P_IDX), 10, 4, 3.5, 'low'); i_Ach = zscore(Filter(i_Ach, 10, 4, bandSensor, 'bandpass'));
    [phi_NE, phi_Ach, spikes, sync, spikes_01, sync_hilbert] = CalcSpikes(i_NE, i_Ach, fs);
%     [RC, PC, LAMBDA, RHO] = SSA(z_pupil, 10);
%     RC1 = RC(:, 1); RC2 = RC(:, 2); RC3 = RC(:, 3);
    
    session_NE{cnt} = z_NE;
    session_Ach{cnt} = z_Ach;
    session_sync{cnt} = sync;
%     [RC, PC, LAMBDA, RHO] = SSA(sync, 20);
%     session_sync_RC1{cnt} = RC(:, 1);
%     session_sync_RC2{cnt} = RC(:, 2);
%     session_sync_RC3{cnt} = RC(:, 3);
    [up, lo] = envelope(sync);

    session_spiketime{cnt} = find(spikes==1);
%     sync = 1-sin(abs(angle(hilbert(i_Ach))-angle(hilbert(i_NE)))/2);
%     m = load(MetaData_files{i});
    mx = load(MetaDataX_files{i});
    
%     if ~isempty(session_Pupil{i})
%         for sf = 1:200
%             rng('shuffle')        % shuffle for sensor signal part
%             session_sensorIDX_shuffle{cnt}(end+1, :) = randperm(length(timestamp));
%         end
%     end
    
    for k = 2:height(mx.MetaDataX)
        if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>3
            if fs==120
                stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Punish_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Tone_Onset(k));
                stTrlOnset_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Trial_Onset(k));
            else
                stOutcome_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Punish_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Tone_Onset(k));
                stTrlOnset_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Trial_Onset(k));
            end
%             session_NE_atOutcome_P{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100);
            session_NE_atOutcome_P{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
%             session_Ach_atOutcome_P{cnt}(end+1, :) = z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100);
            session_Ach_atOutcome_P{cnt}(end+1, :) = z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            session_NA_atOutcome_P{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs)+z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
%             session_NEAch_sync_atOutcome_P{cnt}(end+1, :) = sync(stOutcome_sensor+WINspkcnt_st*fs:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100);
            session_NEAch_sync_atOutcome_P{cnt}(end+1, :) = sync(stOutcome_sensor+WINspkcnt_st*fs:stOutcome_sensor+WINspkcnt_ed*fs);
            session_NEAch_sync_TonePeriod_P{cnt}{end+1, 1} = sync(stTone_sensor:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100); % truncate the trace for TRF pre-OC stim prediction
            session_NEAch_sync_Trial2Outcome_P{cnt}{end+1, 1} = sync(stTrlOnset_sensor:stOutcome_sensor);
%             session_NEAch_syncRC1_TonePeriod_P{cnt}{end+1, 1} = RC(stTone_sensor:stOutcome_sensor, 1);
%             session_NEAch_syncRC2_TonePeriod_P{cnt}{end+1, 1} = RC(stTone_sensor:stOutcome_sensor, 2);
%             session_NEAch_syncRC3_TonePeriod_P{cnt}{end+1, 1} = RC(stTone_sensor:stOutcome_sensor, 3);
            session_NEAch_syncENV_TonePeriod_P{cnt}{end+1, 1} = up(stTone_sensor:stOutcome_sensor);
            session_NE_TonePeriod_P{cnt}{end+1, 1} = z_NE(stTone_sensor:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100); % truncate the trace for TRF pre-OC stim prediction
            session_Ach_TonePeriod_P{cnt}{end+1, 1} = z_Ach(stTone_sensor:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100); % truncate the trace for TRF pre-OC stim prediction
            session_NE_LickFreePeriod_P{cnt}{end+1, 1} = z_NE(stTrlOnset_sensor+WINfree_st*fs:stTone_sensor+WINfree_ed*fs+abs(tmin)/100+abs(tmax)/100);
            session_Ach_LickFreePeriod_P{cnt}{end+1, 1} = z_Ach(stTrlOnset_sensor+WINfree_st*fs:stTone_sensor+WINfree_ed*fs+abs(tmin)/100+abs(tmax)/100);
            session_NEAch_sync_LickFreePeriod_P{cnt}{end+1, 1} = sync(stTrlOnset_sensor+WINfree_st*fs:stTone_sensor+WINfree_ed*fs+abs(tmin)/100+abs(tmax)/100);
            session_spike_TonePeriod_P{cnt}{end+1, 1} = spikes_01(stTone_sensor:stOutcome_sensor);
            session_spike_Trial2Outcome_P{cnt}{end+1, 1} = spikes_01(stTrlOnset_sensor:stOutcome_sensor);
            session_DiffNA_atOutcome_P{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs)-z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            session_sync_atOutcome_NS{cnt}(end+1, :) = sync(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            session_sync_Free_NS{cnt}(end+1, :) = sync(stTone_sensor-10*fs:stTone_sensor+1*fs);
            idx = intersect(find(spikes==1), stOutcome_sensor-5*fs:stOutcome_sensor);
            session_spiketime_P{cnt} = [session_spiketime_P{cnt}; idx];
            session_OutcomeIDX_P{cnt}(end+1, :) = stOutcome_sensor-5*fs:stOutcome_sensor;
            session_spike2P{cnt} = [session_spike2P{cnt}; idx-stOutcome_sensor];
            session_Nspike_P{cnt}(end+1) = length(idx);
%             occupied = [];
%             for sp = 1:length(idx)
%                 occupied = [occupied; [idx-15:idx+15]'];
%             end
%             rng('shuffle');
%             candidate = setdiff(stOutcome_sensor-5*fs:stOutcome_sensor, unique(occupied)); candidate = candidate';
% %             candidate = stOutcome_sensor-4*fs:stOutcome_sensor; candidate = candidate';
%             session_SHFspike_P{cnt} = [session_SHFspike_P{cnt}; candidate(randperm(length(candidate), num_SHF))];
        end
        if mx.MetaDataX.(OPPTWIN{option})(k)==1
            if fs==120
                stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Reward_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Tone_Onset(k));
                stTrlOnset_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Trial_Onset(k));
            else
                stOutcome_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Reward_Onset(k));
                stTone_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Tone_Onset(k));
                stTrlOnset_sensor = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Trial_Onset(k));
            end
%             session_NE_atOutcome_R{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100);
            session_NE_atOutcome_R{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
%             session_Ach_atOutcome_R{cnt}(end+1, :) = z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100);
            session_Ach_atOutcome_R{cnt}(end+1, :) = z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            session_NA_atOutcome_R{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs)+z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
%             session_NEAch_sync_atOutcome_R{cnt}(end+1, :) = sync(stOutcome_sensor+WINspkcnt_st*fs:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100);
            session_NEAch_sync_atOutcome_R{cnt}(end+1, :) = sync(stOutcome_sensor+WINspkcnt_st*fs:stOutcome_sensor+WINspkcnt_ed*fs);
            session_NEAch_sync_TonePeriod_R{cnt}{end+1, 1} = sync(stTone_sensor:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100); % truncate the trace for TRF pre-OC stim prediction);
            session_NEAch_sync_Trial2Outcome_R{cnt}{end+1, 1} = sync(stTrlOnset_sensor:stOutcome_sensor);
%             session_NEAch_syncRC1_TonePeriod_R{cnt}{end+1, 1} = RC(stTone_sensor:stOutcome_sensor, 1);
%             session_NEAch_syncRC2_TonePeriod_R{cnt}{end+1, 1} = RC(stTone_sensor:stOutcome_sensor, 2);
%             session_NEAch_syncRC3_TonePeriod_R{cnt}{end+1, 1} = RC(stTone_sensor:stOutcome_sensor, 3);
            session_NEAch_syncENV_TonePeriod_R{cnt}{end+1, 1} = up(stTone_sensor:stOutcome_sensor);
            session_NE_TonePeriod_R{cnt}{end+1, 1} = z_NE(stTone_sensor:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100); % truncate the trace for TRF pre-OC stim prediction
            session_Ach_TonePeriod_R{cnt}{end+1, 1} = z_Ach(stTone_sensor:stOutcome_sensor+abs(tmin)/100+abs(tmax)/100); % truncate the trace for TRF pre-OC stim prediction
            session_NE_LickFreePeriod_R{cnt}{end+1, 1} = z_NE(stTrlOnset_sensor+WINfree_st*fs:stTone_sensor+WINfree_ed*fs+abs(tmin)/100+abs(tmax)/100);
            session_Ach_LickFreePeriod_R{cnt}{end+1, 1} = z_Ach(stTrlOnset_sensor+WINfree_st*fs:stTone_sensor+WINfree_ed*fs+abs(tmin)/100+abs(tmax)/100);
            session_NEAch_sync_LickFreePeriod_R{cnt}{end+1, 1} = sync(stTrlOnset_sensor+WINfree_st*fs:stTone_sensor+WINfree_ed*fs+abs(tmin)/100+abs(tmax)/100);
            session_spike_TonePeriod_R{cnt}{end+1, 1} = spikes_01(stTone_sensor:stOutcome_sensor);
            session_spike_Trial2Outcome_R{cnt}{end+1, 1} = spikes_01(stTrlOnset_sensor:stOutcome_sensor);
            session_DiffNA_atOutcome_R{cnt}(end+1, :) = z_NE(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs)-z_Ach(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            session_sync_atOutcome_NS{cnt}(end+1, :) = sync(stOutcome_sensor-5*fs:stOutcome_sensor+1*fs);
            session_sync_Free_NS{cnt}(end+1, :) = sync(stTone_sensor-10*fs:stTone_sensor+1*fs);
            idx = intersect(find(spikes==1), stOutcome_sensor-5*fs:stOutcome_sensor);
            session_spiketime_R{cnt} = [session_spiketime_R{cnt}; idx];
            session_OutcomeIDX_R{cnt}(end+1, :) = stOutcome_sensor-5*fs:stOutcome_sensor;
            session_spike2R{cnt} = [session_spike2R{cnt}; idx-stOutcome_sensor];
            session_Nspike_R{cnt}(end+1) = length(idx);
%             occupied = [];
%             for sp = 1:length(idx)
%                 occupied = [occupied; [idx-15:idx+15]'];
%             end
%             rng('shuffle');
%             candidate = setdiff(stOutcome_sensor-5*fs:stOutcome_sensor, unique(occupied)); candidate = candidate';
% %             candidate = stOutcome_sensor-4*fs:stOutcome_sensor; candidate = candidate';
%             session_SHFspike_R{cnt} = [session_SHFspike_R{cnt}; candidate(randperm(length(candidate), num_SHF))];
        end
        
    end
    num0 = size(session_NEAch_sync_atOutcome_P{cnt}, 1);
    num1 = size(session_NEAch_sync_atOutcome_R{cnt}, 1);
%     ldaModel = fitcdiscr([session_NEAch_sync_atOutcome_P{cnt}; session_NEAch_sync_atOutcome_R{cnt}], [zeros(num0, 1); ones(num1, 1)]);
%     transformation_matrix = ldaModel.Coeffs(1, 2).Linear;
%     X_lda = [session_NEAch_sync_atOutcome_P{cnt}; session_NEAch_sync_atOutcome_R{cnt}]*transformation_matrix;
%     session_sync_LDA_P{cnt} = X_lda(1:num0, :);
%     session_sync_LDA_R{cnt} = X_lda(num0+1:end, :);
end
toc
%%
Lumped_NE_atOutcome_P = NaN(length(OFCIDX), size(session_NE_atOutcome_P{1}, 2)); % Punished but withheld >4s NE
Lumped_Ach_atOutcome_P = NaN(length(OFCIDX), size(session_Ach_atOutcome_P{1}, 2)); % Punished but withheld >4s Ach
Lumped_NE_atOutcome_R = NaN(length(OFCIDX), size(session_NE_atOutcome_R{1}, 2)); % Rewarded NE
Lumped_Ach_atOutcome_R = NaN(length(OFCIDX), size(session_Ach_atOutcome_R{1}, 2)); % Rewarded Ach
Lumped_NA_atOutcome_P = NaN(length(OFCIDX), size(session_NEAch_sync_atOutcome_P{1}, 2));
Lumped_NA_atOutcome_R = NaN(length(OFCIDX), size(session_NEAch_sync_atOutcome_R{1}, 2));
for i = 1:length(OFCIDX)
    if isempty(session_NE_atOutcome_P{i}) || isempty(session_Ach_atOutcome_P{i}) || isempty(session_NE_atOutcome_R{i}) || isempty(session_Ach_atOutcome_R{i})...
             || isempty(session_NA_atOutcome_P{i}) || isempty(session_NA_atOutcome_R{i})
        continue;
    end
    Lumped_NE_atOutcome_P(i, :) = nanmean(session_NE_atOutcome_P{i}, 1);
    Lumped_Ach_atOutcome_P(i, :) = nanmean(session_Ach_atOutcome_P{i}, 1);
    Lumped_NE_atOutcome_R(i, :) = nanmean(session_NE_atOutcome_R{i}, 1);
    Lumped_Ach_atOutcome_R(i, :) = nanmean(session_Ach_atOutcome_R{i}, 1);
    Lumped_NA_atOutcome_P(i, :) = nanmean(session_NEAch_sync_atOutcome_P{i}, 1);
    Lumped_NA_atOutcome_R(i, :) = nanmean(session_NEAch_sync_atOutcome_R{i}, 1);
end

subjLumped_NE_atOutcome_P = cell(length(subjIDX), 1); % Punished but withheld >4s NE
subjLumped_Ach_atOutcome_P = cell(length(subjIDX), 1); % Punished but withheld >4s Ach
subjLumped_NE_atOutcome_R = cell(length(subjIDX), 1); % Rewarded NE
subjLumped_Ach_atOutcome_R = cell(length(subjIDX), 1); % Rewarded Ach
subjLumped_NA_atOutcome_P = cell(length(subjIDX), 1);
subjLumped_NA_atOutcome_R = cell(length(subjIDX), 1);
for i = 1:length(subjIDX)
    subjLumped_NE_atOutcome_P{i} = cell2mat(session_NE_atOutcome_P(subjIDX{i}));
    subjLumped_Ach_atOutcome_P{i} = cell2mat(session_Ach_atOutcome_P(subjIDX{i}));
    subjLumped_NE_atOutcome_R{i} = cell2mat(session_NE_atOutcome_R(subjIDX{i}));
    subjLumped_Ach_atOutcome_R{i} = cell2mat(session_Ach_atOutcome_R(subjIDX{i}));
    subjLumped_NA_atOutcome_P{i} = cell2mat(session_NA_atOutcome_P(subjIDX{i}));
    subjLumped_NA_atOutcome_R{i} = cell2mat(session_NA_atOutcome_R(subjIDX{i}));
end

% TT1 = (-2*fs:4*fs)/fs;
% subj_NE_atTone_P = NaN(length(subjIDX), 6*fs+1);
% subj_Ach_atTone_P = NaN(length(subjIDX), 6*fs+1);
%     subj_NE_atTone_R = NaN(length(subjIDX), 6*fs+1);
%     subj_Ach_atTone_R = NaN(length(subjIDX), 6*fs+1);
TT2 = (-5*fs:1*fs)/fs;
subj_NE_atOutcome_P = NaN(length(subjIDX), size(session_NE_atOutcome_P{1}, 2));
subj_Ach_atOutcome_P = NaN(length(subjIDX), size(session_NE_atOutcome_P{1}, 2));
    subj_NE_atOutcome_R = NaN(length(subjIDX), size(session_NE_atOutcome_P{1}, 2));
    subj_Ach_atOutcome_R = NaN(length(subjIDX), size(session_NE_atOutcome_P{1}, 2));
    
for i = 1:length(subjIDX)
%     subj_NE_atTone_P(i, :) = nanmean(Lumped_NE_atTone_P(subjIDX{i}, :), 1);
%     subj_Ach_atTone_P(i, :) = nanmean(Lumped_Ach_atTone_P(subjIDX{i}, :), 1);
%         subj_NE_atTone_R(i, :) = nanmean(Lumped_NE_atTone_R(subjIDX{i}, :), 1);
%         subj_Ach_atTone_R(i, :) = nanmean(Lumped_Ach_atTone_R(subjIDX{i}, :), 1);
    subj_NE_atOutcome_P(i, :) = nanmean(Lumped_NE_atOutcome_P(subjIDX{i}, :), 1);
    subj_Ach_atOutcome_P(i, :) = nanmean(Lumped_Ach_atOutcome_P(subjIDX{i}, :), 1);
        subj_NE_atOutcome_R(i, :) = nanmean(Lumped_NE_atOutcome_R(subjIDX{i}, :), 1);
        subj_Ach_atOutcome_R(i, :) = nanmean(Lumped_Ach_atOutcome_R(subjIDX{i}, :), 1);
end
%%
MEAN_NE_noCNO = [];
MEAN_Ach_noCNO = [];
for i = 1:length(subjIDX)
    MEAN_NE_noCNO(end+1, :) = nanmean([subj_NE_atOutcome_P(i, :); subj_NE_atOutcome_R(i, :)], 1);
    MEAN_Ach_noCNO(end+1, :) = nanmean([subj_Ach_atOutcome_P(i, :); subj_Ach_atOutcome_R(i, :)], 1);
end
%%
MEAN_NE_CNO = [];
MEAN_Ach_CNO = [];
for i = 1:length(subjIDX)
    MEAN_NE_CNO(end+1, :) = nanmean([subj_NE_atOutcome_P(i, :); subj_NE_atOutcome_R(i, :)], 1);
    MEAN_Ach_CNO(end+1, :) = nanmean([subj_Ach_atOutcome_P(i, :); subj_Ach_atOutcome_R(i, :)], 1);
end
%% Plot raw sensor traces before outcome Saline vs CNO
% figure; hold on
subplot(2, 1, 1); 
% ShadedPlot(TT2, mean(MEAN_NE_noCNO, 1), [0 0 1], 1, SEM(MEAN_NE_noCNO), [0.73 0.83 0.96])
ShadedPlot(TT2, mean(MEAN_NE_CNO, 1), [0 0 1], 1, SEM(MEAN_NE_CNO), [0.73 0.83 0.96])
% plot(TT2, mean(MEAN_NE_noCNO, 1), '-b', 'LineWidth', 1)
plot(TT2, mean(MEAN_NE_CNO, 1), '-c', 'LineWidth', 1)
ylabel('Z-score');
vline(0, '-k')

% figure; hold on
subplot(2, 1, 2); 
% ShadedPlot(TT2, mean(MEAN_Ach_noCNO, 1), [0 0 1], 1, SEM(MEAN_Ach_noCNO), [0.9 0.8 0.7])
ShadedPlot(TT2, mean(MEAN_Ach_CNO, 1), [0 0 1], 1, SEM(MEAN_Ach_CNO), [0.9 0.8 0.7])
% plot(TT2, mean(MEAN_Ach_noCNO, 1), '-r', 'LineWidth', 1)
plot(TT2, mean(MEAN_Ach_CNO, 1), '-m', 'LineWidth', 1)
ylabel('Z-score');
xlabel('Time after Outcome (s)')
vline(0, '-k')

%%


%%
function FR = CalcFR(SPK, ONSET, EDGE, STEPSZ)
FR = NaN(1, length(EDGE)-1);
for i = 1:length(EDGE)-1
    FR(i) = length(find(~isnan(SPK(ONSET+EDGE(i)*10:ONSET+EDGE(i+1)*10))));
end
FR = FR/STEPSZ;
end