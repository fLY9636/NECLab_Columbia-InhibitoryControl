% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
    MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);

Performance_noCNO_DBh
Performance_CNO_DBh
SIMPerformance_noCNO_DBh

% Performance_noCNO_Th
% Performance_CNO_Th
% SIMPerformance_noCNO_Th
% 
% Performance_noCNO_ChAT
% Performance_CNO_ChAT
% SIMPerformance_noCNO_ChAT

% Performance_noCNO_OFC
% SIMPerformance_noCNO_OFC

% Performance_noCNO_WTPPC
% SIMPerformance_noCNO_WTPPC

% Performance_noCNO_WTS1
% SIMPerformance_noCNO_WTS1

% Performance_noCNO_WTThm
% SIMPerformance_noCNO_WTThm

SuccessRate = Performance_CNO_DBh.Outcome1(:, 1);
CL_SuccessRate = SIMPerformance_noCNO_DBh;
%% Calculate the AUROC of signals (NE, ACh, encoder value, pupil, HPF pupil)
% session_Pupil_AUROC = NaN(len, 1);
session_NE_AUROC = NaN(len, 1);
session_Ach_AUROC = NaN(len, 1);
session_sync_AUROC = NaN(len, 1);
WINst_pupil = 2*10+1;
WINed_pupil = 5*10;
WINst_sensor = 2*10+1;
WINed_sensor = 4*10;
cnt = 0;
for i = OFCIDX(:).'
    cnt = cnt+1;
    dF_NE = nanmean(session_NE_atOutcome_P{cnt}(:, WINst_sensor:WINed_sensor), 2);
    dS_NE = nanmean(session_NE_atOutcome_R{cnt}(:, WINst_sensor:WINed_sensor), 2);
    dF_Ach = nanmean(session_Ach_atOutcome_P{cnt}(:, WINst_sensor:WINed_sensor), 2);
    dS_Ach = nanmean(session_Ach_atOutcome_R{cnt}(:, WINst_sensor:WINed_sensor), 2);
    dF_NA = nanmean(session_NEAch_sync_atOutcome_P{cnt}(:, WINst_sensor:WINed_sensor), 2);
    dS_NA = nanmean(session_NEAch_sync_atOutcome_R{cnt}(:, WINst_sensor:WINed_sensor), 2);
    session_NE_AUROC(cnt) = abs(calculate_roc_auc_manual(dF_NE, dS_NE));
    session_Ach_AUROC(cnt) = abs(calculate_roc_auc_manual(dF_Ach, dS_Ach));
    session_sync_AUROC(cnt) = abs(calculate_roc_auc_manual(dF_NA, dS_NA));
%     if isempty(session_HiPupil_atOutcome_P{cnt}) || isempty(session_HiPupil_atOutcome_R{cnt}) 
%         continue;
%     end
%     dF = nanmean(session_HiPupil_atOutcome_P{cnt}(:, WINst_pupil:WINed_pupil), 2);
%     dS = nanmean(session_HiPupil_atOutcome_R{cnt}(:, WINst_pupil:WINed_pupil), 2);
%     session_Pupil_AUROC(cnt) = abs(calculate_roc_auc_manual(dF, dS));
end

%%
session_NA_spikeratio = NaN(len, 1); % Punished but withheld >4s NE
session_sync_dprime = NaN(len, 1);
session_NA_spikediff = NaN(len, 1); % Punished but withheld >4s NE
session_PC_P = NaN(len, 1);
session_PC_R = NaN(len, 1);
session_NA_syncdiff = NaN(len, 1);

tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 3; 
bandSensor = [0.4 0.8];
WINspkcnt_st = -3;
WINspkcnt_ed = 0;

cnt = 0;
% figure; hold on
for i = OFCIDX(:).'
    cnt = cnt+1; 
    S2P_IDX = e1.MetaData.Pupil_in_sensor{i};
    z_NE = zscore(Filter(e1.MetaData.NE_470{i}, 120, 2, [0.1 10], 'bandpass'));
    z_Ach = zscore(Filter(e1.MetaData.Ach_470{i}, 120, 2, [0.1 10], 'bandpass'));

    i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, bandSensor, 'bandpass');
    i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, bandSensor, 'bandpass');
    [phi_NE, phi_Ach, spikes, sync, spikes_01, sync_hilbert] = CalcSpikes(i_NE, i_Ach, 120);
    
    sessiontrial_spikecnt_P = [];
    sessiontrial_spikecnt_R = [];
%     session_NEAch_sync_atOutcome_P = [];
%     session_NEAch_sync_atOutcome_R = [];
    sessiontrial_ave_sync_P = [];
    sessiontrial_ave_sync_R = [];
    sessiontrial_PC_P = [];
    sessiontrial_PC_R = [];
    mx = load(MetaDataX_files{i});
 
    for k = 1:height(mx.MetaDataX)
        if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>3 % (CHANGE starting time back to '3' for RUNNING 'extractPupil.m')
            
            stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Punish_Onset(k));
            
            stOutcome_spike = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Punish_Onset(k));
            sessiontrial_spikecnt_P(end+1, 1) = length(find(~isnan(spikes(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120))))/(WINspkcnt_ed-WINspkcnt_st);
%             session_NEAch_sync_atOutcome_P{cnt}(end+1, :) = sync(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120);
            sessiontrial_ave_sync_P(end+1, :) = nanmean(sync(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120));
            sessiontrial_PC_P = [sessiontrial_PC_P; calculate_phase_coherence(phi_NE(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120), ...
                phi_Ach(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120), 120, 0.2)];
        end
        if mx.MetaDataX.(OPPTWIN{option})(k)==1
            stOutcome_sensor = CrossSampling(e1.MetaData.Sensor_time{i}, mx.MetaDataX.Reward_Onset(k));
            
            stOutcome_spike = CrossSampling(e1.MetaData.Pupil_time{i}, mx.MetaDataX.Reward_Onset(k));
            sessiontrial_spikecnt_R(end+1, 1) = length(find(~isnan(spikes(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120))))/(WINspkcnt_ed-WINspkcnt_st);
%             session_NEAch_sync_atOutcome_R{cnt}(end+1, :) = sync(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120);
            sessiontrial_ave_sync_R(end+1, :) = nanmean(sync(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120));
            sessiontrial_PC_R = [sessiontrial_PC_R; calculate_phase_coherence(phi_NE(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120), ...
                phi_Ach(stOutcome_sensor+WINspkcnt_st*120:stOutcome_sensor+WINspkcnt_ed*120), 120, 0.2)];
        end
    end
    session_NA_spikeratio(cnt) = mean(sessiontrial_spikecnt_P)/mean(sessiontrial_spikecnt_R);
    session_NA_spikediff(cnt) = mean(sessiontrial_spikecnt_P)-mean(sessiontrial_spikecnt_R);
    session_NA_syncdiff(cnt) = mean(sessiontrial_ave_sync_P)-mean(sessiontrial_ave_sync_R);
    if isempty(sessiontrial_ave_sync_R) || isempty(sessiontrial_ave_sync_P)
        session_sync_dprime(cnt) = NaN;
    else
        session_sync_dprime(cnt) = norm(nanmean(sessiontrial_ave_sync_P, 1)-nanmean(sessiontrial_ave_sync_R, 1));
    end
    session_PC_P(cnt) = nanmean(sessiontrial_PC_P);
    session_PC_R(cnt) = nanmean(sessiontrial_PC_R);
%     session_PC(cnt) = calculate_phase_coherence(sessiontrial_phi_NE, sessiontrial_phi_Ach, );
%     subplot(7, 7, cnt); hold on
%     histogram(sessiontrial_ave_sync_P, 10); 
%     histogram(sessiontrial_ave_sync_R, 10); 
%     title([num2str(nanmean(sessiontrial_ave_sync_P)-nanmean(sessiontrial_ave_sync_R))])
end
toc
%%
Ratio_Saline = session_NA_spikeratio;
dprime_Saline = session_sync_dprime;
Diff_Saline = session_NA_spikediff;
syncDiff_Saline = session_NA_syncdiff;
% PC_Saline = session_PC;
%%
Ratio_CNO = session_NA_spikeratio;
dprime_CNO = session_sync_dprime;
Diff_CNO = session_NA_spikediff;
syncDiff_CNO = session_NA_syncdiff;
% PC_CNO = session_PC;
%%
%%%%% session baseline norm
norm_SuccessRate = SuccessRate;
for i = 1:length(subjIDX)
    norm_SuccessRate(subjIDX{i}) = CalcPercent(norm_SuccessRate(subjIDX{i}), nanmean(CL_SuccessRate{i}));
end
SuccessRate = norm_SuccessRate;
%%%%% session baseline diff
% norm_SuccessRate = SuccessRate;
% for i = 1:length(subjIDX)
%     norm_SuccessRate(subjIDX{i}) = norm_SuccessRate(subjIDX{i})-mean(CL_SuccessRate{i});
% end
% SuccessRate = norm_SuccessRate*100;
%%%%% subj MIN-MAX norm
% norm_SuccessRate = SuccessRate;
% for i = 1:length(subjIDX)
%     norm_SuccessRate(subjIDX{i}) = UnitNormalization(norm_SuccessRate(subjIDX{i}));
% end
% SuccessRate = norm_SuccessRate*100;
%%
SuccessRate_Saline = SuccessRate;
%%
SuccessRate_CNO = SuccessRate;
%% Use polyfit[1] (use this for DBH syncDiff_Saline , syncDiff_CNO vs Performance)
% Assuming you have two vectors X and Y
% X = [SuccessRate_Saline; SuccessRate_CNO];  % Replace this with your X data
% Y = [Ratio_Saline; Ratio_CNO];  % Example: Y = 2*X + 1 with some random noise
X = SuccessRate_Saline;  % Replace this with your X data
Y= syncDiff_Saline;  % Example: Y = 2*X + 1 with some random noise
% X = SuccessRate_CNO;  % Replace this with your X data
% Y = Diff_CNO;  % Example: Y = 2*X + 1 with some random noise
% Create the X-Y scatter plot
valid_indices = isfinite(X) & isfinite(Y);
X_valid = X(valid_indices);
Y_valid = Y(valid_indices);
subplot(2, 1, 2); hold on
scatter(X, Y, 12);  % 'b' specifies the color as blue
% scatter(X, Y, 'c', 'filled');  % 'b' specifies the color as blue
% Fit the linear line usng polyfit
coefficients = polyfit(X_valid, Y_valid, 1);
% Get the slope and y-intercept of the fitted line
slope = coefficients(1);
intercept = coefficients(2);
% Create the fitted line using polyval
fitted_line = slope*X_valid+intercept;
% Plot the fitted line
hold on;
plot(X_valid, fitted_line, 'r', 'LineWidth', 2);  % 'r' specifies the color as red
% Add labels and title
xlabel('Success rate %');
ylabel('Switching rate difference \newline [failed -- successful]');
title(['Linear fit [' num2str(WINspkcnt_st) 's, ' num2str(WINspkcnt_ed) 's] OC']);
% xlim([0 100])
% Calculate the coefficient of determination (R-squared)
Y_mean = mean(Y_valid);
SS_total = sum((Y_valid-Y_mean).^2);
SS_residual = sum((Y_valid-fitted_line).^2);
R_squared = 1-SS_residual/SS_total;
[p, ~, mu] = polyfit(X_valid, Y_valid, 1);
% Display the statistics
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);
disp(['R-squared: ', num2str(R_squared)]);
disp(['p Value: ', num2str(p)]);
%% Use fitlm
figure
% Assuming you have two vectors X and Y
% X = [SuccessRate_Saline; SuccessRate_CNO];  % Replace this with your X data
% Y = [Ratio_Saline; Ratio_CNO];  % Example: Y = 2*X + 1 with some random noise
% X = SuccessRate_Saline;  % Replace this with your X data
% Y = Diff_Saline;  % Example: Y = 2*X + 1 with some random noise
X = SuccessRate_CNO;  % Replace this with your X data
Y = -syncDiff_CNO;  % Example: Y = 2*X + 1 with some random noise
% Create the X-Y scatter plot
valid_indices = isfinite(X) & isfinite(Y);
X_valid = X(valid_indices);
Y_valid = Y(valid_indices);
% subplot(2, 1, 2); hold on
scatter(X, Y, 12);  % 'b' specifies the color as blue
% scatter(X, Y, 'c', 'filled');  % 'b' specifies the color as blue
% Fit the linear line usng polyfit
lm = fitlm(X_valid, Y_valid);
% Get the slope and y-intercept of the fitted line
slope = lm.Coefficients.Estimate(2);
intercept = lm.Coefficients.Estimate(1);
% Create the fitted line using polyval
fitted_values = predict(lm, X_valid);
% Plot the fitted line
hold on;
plot(X_valid, fitted_values, 'r', 'LineWidth', 2);  % 'r' specifies the color as red
% Add labels and title
xlabel('Normalized success rate %');
% xlabel('AUROC NE');
% xlabel('Switching rate difference \newline [failed -- successful]');
% ylabel('Switching rate difference \newline [failed -- successful]');
ylabel('Phase synchrony difference \newline [successful -- failed]');
% ylabel('AUROC HPF pupil');
title(['NE']);
% xlim([0 100])
% Calculate the coefficient of determination (R-squared)
% Y_mean = mean(Y_valid);
% SS_total = sum((Y_valid-Y_mean).^2);
% SS_residual = sum((Y_valid-fitted_line).^2);
R_squared = lm.Rsquared.Ordinary;
p_value = lm.Coefficients.pValue(2);
% Display the statistics
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);
disp(['R-squared: ', num2str(R_squared)]);
disp(['p Value: ', num2str(p_value)]);
%% Plot the distribution (or raw trace, for pupil) of a specific session based on different AUROC level
X = session_Pupil_AUROC;
X_P = session_Pupil_atOutcome_P;
X_R = session_Pupil_atOutcome_R;
figure; hold on
IDX_HiAUROC = find(X>=0.7);
X_P_temp = cell2mat(X_P(IDX_HiAUROC));
X_R_temp = cell2mat(X_R(IDX_HiAUROC));
Distr_P = nanmean(X_P_temp(:, WINst_pupil:WINed_pupil), 2);
Distr_R = nanmean(X_R_temp(:, WINst_pupil:WINed_pupil), 2);
histogram(Distr_P)
histogram(Distr_R)
figure; hold on
ShadedPlot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_P_temp(:, WINst_pupil:WINed_pupil), 1), 'k', 1, SEM(X_P_temp(:, WINst_pupil:WINed_pupil)), [0.8 0.8 0.8])
ShadedPlot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_R_temp(:, WINst_pupil:WINed_pupil), 1), 'k', 1, SEM(X_R_temp(:, WINst_pupil:WINed_pupil)), [0.8 0.8 0.8])
plot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_P_temp(:, WINst_pupil:WINed_pupil), 1), '--k')
plot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_R_temp(:, WINst_pupil:WINed_pupil), 1), '-k')

figure; hold on
IDX_LoAUROC = find(X<=0.35);
X_P_temp = cell2mat(X_P(IDX_LoAUROC));
X_R_temp = cell2mat(X_R(IDX_LoAUROC));
Distr_P = nanmean(X_P_temp(:, WINst_pupil:WINed_pupil), 2);
Distr_R = nanmean(X_R_temp(:, WINst_pupil:WINed_pupil), 2);
histogram(Distr_P)
histogram(Distr_R)
figure; hold on
ShadedPlot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_P_temp(:, WINst_pupil:WINed_pupil), 1), 'k', 1, SEM(X_P_temp(:, WINst_pupil:WINed_pupil)), [0.8 0.8 0.8])
ShadedPlot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_R_temp(:, WINst_pupil:WINed_pupil), 1), 'k', 1, SEM(X_R_temp(:, WINst_pupil:WINed_pupil)), [0.8 0.8 0.8])
plot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_P_temp(:, WINst_pupil:WINed_pupil), 1), '--k')
plot((WINst_pupil:WINed_pupil)/10-5, nanmean(X_R_temp(:, WINst_pupil:WINed_pupil), 1), '-k')





%%
X1 = Performance_noCNO_Th.Outcome1(:, 1)*100;  % Replace this with your X data
Y1= Ratio_Saline;  % Example: Y = 2*X + 1 with some random noise
valid_indices = isfinite(X1) & isfinite(Y1);
X_valid = X1(valid_indices);
Y_valid = Y1(valid_indices);
coefficients = polyfit(X_valid, Y_valid, 1);
slope = coefficients(1);
intercept = coefficients(2);
fitted_line = slope*X_valid+intercept;
Y_mean = mean(Y_valid);
SS_total = sum((Y_valid-Y_mean).^2);
SS_residual = sum((Y_valid-fitted_line).^2);
R_squared = 1-SS_residual/SS_total;
[p, ~, mu] = polyfit(X_valid, Y_valid, 1);
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);
disp(['R-squared: ', num2str(R_squared)]);
disp(['p Value: ', num2str(p)]);










