%% Create the nested data structure for NE, ACh or encoder value
TRFcmp_ACh = struct;
TRFcmp_ACh.cv_corr_P = Lumped_cv_corr_P;
TRFcmp_ACh.cv_corr_R = Lumped_cv_corr_R;
TRFcmp_ACh.cv_mse_P = Lumped_cv_mse_P;
TRFcmp_ACh.cv_mse_R = Lumped_cv_mse_R;
TRFcmp_ACh.pred_corr_P = Lumped_pred_corr_P;
TRFcmp_ACh.pred_corr_R = Lumped_pred_corr_R;
TRFcmp_ACh.pred_stim_P = Lumped_pred_stim_P;
TRFcmp_ACh.pred_stim_R = Lumped_pred_stim_R;
TRFcmp_ACh.w_coef_P = Lumped_w_coeff_P;
TRFcmp_ACh.w_coef_R = Lumped_w_coeff_R;
TRFcmp_ACh.wtrans_coef_P = Lumped_wtrans_coeff_P;
TRFcmp_ACh.wtrans_coef_R = Lumped_wtrans_coeff_R;
TRFcmp_ACh.model_bwd = model_bwd;
TRFcmp_ACh.model_fwd = model_fwd;
TRFcmp_ACh.model_trans = model_trans;
%%
TRFcmp_NE = struct;
TRFcmp_NE.cv_corr_P = Lumped_cv_corr_P;
TRFcmp_NE.cv_corr_R = Lumped_cv_corr_R;
TRFcmp_NE.cv_mse_P = Lumped_cv_mse_P;
TRFcmp_NE.cv_mse_R = Lumped_cv_mse_R;
TRFcmp_NE.pred_corr_P = Lumped_pred_corr_P;
TRFcmp_NE.pred_corr_R = Lumped_pred_corr_R;
TRFcmp_NE.pred_stim_P = Lumped_pred_stim_P;
TRFcmp_NE.pred_stim_R = Lumped_pred_stim_R;
TRFcmp_NE.w_coef_P = Lumped_w_coeff_P;
TRFcmp_NE.w_coef_R = Lumped_w_coeff_R;
TRFcmp_NE.wtrans_coef_P = Lumped_wtrans_coeff_P;
TRFcmp_NE.wtrans_coef_R = Lumped_wtrans_coeff_R;
TRFcmp_NE.model_bwd = model_bwd;
TRFcmp_NE.model_fwd = model_fwd;
TRFcmp_NE.model_trans = model_trans;
%%
TRFcmp_sync = struct;
TRFcmp_sync.cv_corr_P = Lumped_cv_corr_P;
TRFcmp_sync.cv_corr_R = Lumped_cv_corr_R;
TRFcmp_sync.cv_mse_P = Lumped_cv_mse_P;
TRFcmp_sync.cv_mse_R = Lumped_cv_mse_R;
TRFcmp_sync.pred_corr_P = Lumped_pred_corr_P;
TRFcmp_sync.pred_corr_R = Lumped_pred_corr_R;
TRFcmp_sync.pred_stim_P = Lumped_pred_stim_P;
TRFcmp_sync.pred_stim_R = Lumped_pred_stim_R;
TRFcmp_sync.w_coef_P = Lumped_w_coeff_P;
TRFcmp_sync.w_coef_R = Lumped_w_coeff_R;
TRFcmp_sync.wtrans_coef_P = Lumped_wtrans_coeff_P;
TRFcmp_sync.wtrans_coef_R = Lumped_wtrans_coeff_R;
TRFcmp_sync.model_bwd = model_bwd;
TRFcmp_sync.model_fwd = model_fwd;
TRFcmp_sync.model_trans = model_trans;
%% Plot the reconstructed traces
X = TRFcmp_NE;
MEAN_pred_stim_P = NaN(count, 30);
MEAN_pred_stim_R = NaN(count, 30);
for i = 1:count
    MEAN_pred_stim_P(i, :) = nanmean(X.pred_stim_P{i}, 1);
    MEAN_pred_stim_R(i, :) = nanmean(X.pred_stim_R{i}, 1);
end
figure; hold on
ShadedPlot((1:30)/10-3, nanmean(MEAN_pred_stim_P, 1), 'b', 1, SEM(MEAN_pred_stim_P), [0.8 0.8 0.8])
ShadedPlot((1:30)/10-3, nanmean(MEAN_pred_stim_R, 1), 'b', 1, SEM(MEAN_pred_stim_R), [0.8 0.8 0.8])
plot((1:30)/10-3, nanmean(MEAN_pred_stim_P, 1), '--b', 'LineWidth', 1);
plot((1:30)/10-3, nanmean(MEAN_pred_stim_R, 1), '-b', 'LineWidth', 1);
xlabel('Time after outcome (s)')
ylabel('Predicted phase synchrony')
%% Plot the mean predicted encoder value [-3s, -1s] before outcome
figure; hold on
bar(1, nanmean(nanmean(MEAN_pred_stim_P(:, 1:20), 2)))
bar(2, nanmean(nanmean(MEAN_pred_stim_R(:, 1:20), 2)))
for i = 1:count
    plot([1 2], [nanmean(MEAN_pred_stim_P(i, 1:20), 2) nanmean(MEAN_pred_stim_R(i, 1:20), 2)], '-ok')
end
errorbar([1 2], [nanmean(nanmean(MEAN_pred_stim_P(:, 1:20), 2)) nanmean(nanmean(MEAN_pred_stim_R(:, 1:20), 2))], ...
    [SEM(nanmean(MEAN_pred_stim_P(:, 1:20), 2)) SEM(nanmean(MEAN_pred_stim_R(:, 1:20), 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean predicted phase synchrony')
%% Plot the correlation score of test set
MEAN_pred_corr_P = NaN(count, 1);
MEAN_pred_corr_R = NaN(count, 1);
for i = 1:count
    MEAN_pred_corr_P(i) = nanmean(X.pred_corr_P{i});
    MEAN_pred_corr_R(i) = nanmean(X.pred_corr_R{i});
end
figure; hold on
bar(1, nanmean(MEAN_pred_corr_P))
bar(2, nanmean(MEAN_pred_corr_R))
for i = 1:count
    plot([1:2], [MEAN_pred_corr_P(i) MEAN_pred_corr_R(i)], '-ok')
end
errorbar([1 2], [nanmean(MEAN_pred_corr_P) nanmean(MEAN_pred_corr_R)], ...
    [SEM(MEAN_pred_corr_P) SEM(MEAN_pred_corr_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Testing performance [r]')
%% Plot the filter weights across time lags
X = TRFcmp_sync;
MEAN_w_P = [];
MEAN_w_R = [];
ZMEAN_w_P = [];
ZMEAN_w_R = [];
sessionMEAN_w_P = NaN(count, length(X.model_fwd.t));
sessionMEAN_w_R = NaN(count, length(X.model_fwd.t));
sessionZMEAN_w_P = NaN(count, length(X.model_fwd.t));
sessionZMEAN_w_R = NaN(count, length(X.model_fwd.t));
for i = 1:count
    Pcell = X.w_coef_P{i};
    Rcell = X.w_coef_R{i};
    MEAN_w_P = [MEAN_w_P; Pcell];
    MEAN_w_R = [MEAN_w_R; Rcell];
%     ZMEAN_w_P = [ZMEAN_w_P; zscore(Pcell, [], 2)];
    ZMEAN_w_P = [ZMEAN_w_P; Pcell-nanmean(Pcell, 2)];
%     ZMEAN_w_R = [ZMEAN_w_R; zscore(Rcell, [], 2)];
    ZMEAN_w_R = [ZMEAN_w_R; Rcell-nanmean(Rcell, 2)];
    sessionMEAN_w_P(i, :) = nanmean(Pcell, 1);
    sessionMEAN_w_R(i, :) = nanmean(Rcell, 1);
%     sessionZMEAN_w_P(i, :) = nanmean(zscore(Pcell, [], 2), 1);
    sessionZMEAN_w_P(i, :) = nanmean(Pcell-nanmean(Pcell, 2), 1);
%     sessionZMEAN_w_R(i, :) = nanmean(zscore(Rcell, [], 2), 1);
    sessionZMEAN_w_R(i, :) = nanmean(Rcell-nanmean(Rcell, 2), 1);
end


figure; hold on
ShadedPlot(X.model_fwd.t/1000, nanmean(MEAN_w_P, 1), 'k', 1, SEM(MEAN_w_P), [0.8 0.8 0.8])
ShadedPlot(X.model_fwd.t/1000, nanmean(MEAN_w_R, 1), 'k', 1, SEM(MEAN_w_R), [0.8 0.8 0.8])
plot(X.model_fwd.t/1000, nanmean(MEAN_w_P, 1), '--k', 'LineWidth', 1);
plot(X.model_fwd.t/1000, nanmean(MEAN_w_R, 1), '-k', 'LineWidth', 1);
vline(0, ':k')
xlabel('Time lag (s) [sensor - pupil]')
ylabel('a.u.')

figure; hold on
ShadedPlot(X.model_fwd.t/1000, nanmean(ZMEAN_w_P, 1), 'k', 1, SEM(ZMEAN_w_P), [0.8 0.8 0.8])
ShadedPlot(X.model_fwd.t/1000, nanmean(ZMEAN_w_R, 1), 'k', 1, SEM(ZMEAN_w_R), [0.8 0.8 0.8])
plot(X.model_fwd.t/1000, nanmean(ZMEAN_w_P, 1), '--k', 'LineWidth', 1);
plot(X.model_fwd.t/1000, nanmean(ZMEAN_w_R, 1), '-k', 'LineWidth', 1);
vline(0, ':k')
xlabel('Time lag (s) [sensor - pupil]')
ylabel('Z a.u.')

W_P = X.w_coef_P;
W_R = X.w_coef_R;
pos_w_IDXst = CrossSampling(X.model_fwd.t, 0);
pos_w_IDXed = CrossSampling(X.model_fwd.t, 1000); % 600 for ave; 800 for area
signfcnt_T_P = 0.1*ones(1, length(X.model_fwd.t));
signfcnt_T_R = -0.1*ones(1, length(X.model_fwd.t));
signfcnt_T = zeros(1, length(X.model_fwd.t));
for tt = 1:length(signfcnt_T)
    [h, p, t] = qw_statPairedTest(sessionMEAN_w_P(:, tt), 0);
    if p <= 0.05
        signfcnt_T_P(tt) = 1;
%         display(p)
    end
    [h, p, t] = qw_statPairedTest(sessionMEAN_w_R(:, tt), 0);
    if p <= 0.05
        signfcnt_T_R(tt) = -1;
%         display(p)
    end
    [h, p, t] = qw_statPairedTest(sessionZMEAN_w_P(:, tt)-sessionZMEAN_w_R(:, tt), 0);
    if p <= 0.05
        signfcnt_T(tt) = 1.3;
%                 display(p)
    end
end
figure; hold on
scatter(X.model_fwd.t/1000, signfcnt_T_P);
scatter(X.model_fwd.t/1000, signfcnt_T_R);
scatter(X.model_fwd.t/1000, signfcnt_T);

Pvalid = find(~isnan(sessionMEAN_w_P(:, 1)));
Rvalid = find(~isnan(sessionMEAN_w_R(:, 1)));
% Lumped_ZMEAN_w = [sessionZMEAN_w_P(Pvalid, :); sessionZMEAN_w_R(Rvalid, :)];
% [p, tbl] = anova2(Lumped_ZMEAN_w, 1);
data_P = cell(1, length(signfcnt_T)-3);
data_R = cell(1, length(signfcnt_T)-3);
group_w_Pvalid = sessionMEAN_w_P(Pvalid, :);
group_w_Rvalid = sessionMEAN_w_R(Rvalid, :);
cnt1 = 0;
for tt = 4:length(signfcnt_T)
    cnt1 = cnt1+1;
    data_P{cnt1} = group_w_Pvalid(:, tt)';
    data_R{cnt1} = group_w_Rvalid(:, tt)';
end
% performANOVA(data_P)
% performANOVA(data_R)

[h, p, t] = qw_statPairedTest(nanmean(sessionZMEAN_w_P(:, 4:end), 2), nanmean(sessionZMEAN_w_R(:, 4:end), 2))

%% Plot the filter weights across time lags (transformed from the backward model)
MEAN_w_P = [];
MEAN_w_R = [];
sessionMEAN_wtrans_P = NaN(count, length(model_fwd.t));
sessionMEAN_wtrans_R = NaN(count, length(model_fwd.t));
for i = 1:count
    MEAN_w_P = [MEAN_w_P; X.wtrans_coef_P{i}];
    MEAN_w_R = [MEAN_w_R; X.wtrans_coef_R{i}];
    sessionMEAN_wtrans_P(i, :) = nanmean(X.wtrans_coef_P{i}, 1);
    sessionMEAN_wtrans_R(i, :) = nanmean(X.wtrans_coef_R{i}, 1);
end
figure; hold on
ShadedPlot(model_fwd.t/1000, nanmean(MEAN_w_P, 1), 'k', 1, SEM(MEAN_w_P), [0.8 0.8 0.8])
ShadedPlot(model_fwd.t/1000, nanmean(MEAN_w_R, 1), 'k', 1, SEM(MEAN_w_R), [0.8 0.8 0.8])
plot(model_fwd.t/1000, nanmean(MEAN_w_P, 1), '--k', 'LineWidth', 1);
plot(model_fwd.t/1000, nanmean(MEAN_w_R, 1), '-k', 'LineWidth', 1);
vline(0, ':k')
xlabel('Time lag (s) [NE - HPF pupil]')
ylabel('a.u.')
%%
session_sync_ave_pred_corr = NaN(96, 1);
session_NE_ave_pred_corr = NaN(96, 1);
session_Ach_ave_pred_corr = NaN(96, 1);
for i = 1:count
    session_sync_ave_pred_corr(i) = (nanmean(TRFcmp_sync.pred_corr_P{i})+nanmean(TRFcmp_sync.pred_corr_R{i}))/2;
    session_NE_ave_pred_corr(i) = (nanmean(TRFcmp_NE.pred_corr_P{i})+nanmean(TRFcmp_NE.pred_corr_R{i}))/2;
    session_Ach_ave_pred_corr(i) = (nanmean(TRFcmp_Ach.pred_corr_P{i})+nanmean(TRFcmp_Ach.pred_corr_R{i}))/2;
end
figure; hold on
for i = 1:count
    plot([1 2 3], [session_NE_ave_pred_corr(i) session_Ach_ave_pred_corr(i) session_sync_ave_pred_corr(i)], '-o', 'Color', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8])
end
plot([1 2 3], [nanmean(session_NE_ave_pred_corr) nanmean(session_Ach_ave_pred_corr) nanmean(session_sync_ave_pred_corr)], '-om', 'LineWidth', 1.5)
errorbar([1 2 3], [nanmean(session_NE_ave_pred_corr) nanmean(session_Ach_ave_pred_corr) nanmean(session_sync_ave_pred_corr)], ...
    [SEM(session_NE_ave_pred_corr) SEM(session_Ach_ave_pred_corr) SEM(session_sync_ave_pred_corr)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

%%
SuccessRate = Performance_noCNO_WT.Outcome1(:, 1);
CL_SuccessRate = SIMPerformance_noCNO_WT;
norm_SuccessRate = SuccessRate;
for i = 1:length(subjIDX)
    norm_SuccessRate(subjIDX{i}) = CalcPercent(norm_SuccessRate(subjIDX{i}), nanmean(CL_SuccessRate{i}));
end
SuccessRate = norm_SuccessRate;
figure; hold on
cnt = 0;
for i = 1:len
    if isempty(session_Pupil{i})
        continue
    end
    cnt = cnt+1;
    subplot(8, 12, cnt); hold on
    title([num2str(round(SuccessRate(i))) ' %'])
    plot(TRFcmp_sync.model_fwd.t/1000, nanmean(TRFcmp_sync.w_coef_P{cnt}, 1), 'LineWidth', 1.5)
end
%%
figure; hold on
cnt = 0;
for i = 1:len
    if isempty(session_Pupil{i})
        continue
    end
    cnt = cnt+1;
    subplot(8, 12, cnt); hold on
    title([num2str(round(SuccessRate(i))) ' %'])
    plot(TRFcmp_sync.model_fwd.t/1000, nanmean(TRFcmp_sync.w_coef_R{cnt}, 1), 'LineWidth', 1.5)
end

%%
pos_w_IDXst = CrossSampling(X.model_fwd.t, 0);
pos_w_IDXed = CrossSampling(X.model_fwd.t, 600); % 600 for ave; 800 for area
Data = sessionMEAN_wtrans_R(:, pos_w_IDXst:pos_w_IDXed);
[coeff, score, latent, tsquared, explained, mu] = pca(Data);
% GMModel = fitgmdist(Data, 2);
% figure; hold on
scatter3(score(:, 1), score(:, 2), score(:, 3))
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')

%%
X = TRFcmp_Ach;
trialMEAN_w = cell(count, 1);
OC = cell(count, 1);
W_P = X.wtrans_coef_P;
W_R = X.wtrans_coef_R;
pos_w_IDXst = CrossSampling(X.model_fwd.t, 0);
pos_w_IDXed = CrossSampling(X.model_fwd.t, 600); % 600 for ave; 800 for area
for i = 1:count
    trialMEAN_w{i} = [trialMEAN_w{i}; nanmean(W_P{i}(:, pos_w_IDXst:pos_w_IDXed), 2); nanmean(W_R{i}(:, pos_w_IDXst:pos_w_IDXed), 2)];
    OC{i} = [OC{i}; zeros(length(nanmean(W_P{i}, 2)), 1); ones(length(nanmean(W_R{i}, 2)), 1)];
    valid_idx = find(~isnan(trialMEAN_w{i}));
%     trialMEAN_w{i} = zscore(trialMEAN_w{i}(valid_idx));
    trialMEAN_w{i} = trialMEAN_w{i}(valid_idx);
    OC{i} = OC{i}(valid_idx);
end
%%
num_bins = 20;
SuccessRate_Binned = NaN(count, num_bins);
QT_Binned = NaN(count, num_bins);
Slope_set = NaN(count, 1);
for i = 1:count
    [MeanINBin, SEMINBin, bin_ctrs, bin_IDXs] = BinnedQTMean1d(UnitNormalization(trialMEAN_w{i}), OC{i}, num_bins);
    for j = 1:num_bins
        num_trials = length(bin_IDXs{j});
        if num_trials==0
            continue;
        end
        num_success = length(find(OC{i}(bin_IDXs{j})==1));
        SuccessRate_Binned(i, j) = num_success/num_trials;
    end
    QT_Binned(i, :) = bin_ctrs-1/num_bins/2;

    XX = bin_ctrs-1/num_bins/2;
    YY = SuccessRate_Binned(i, :);
    lm = fitlm(XX, YY);
    Slope_set(i) = lm.Coefficients.Estimate(2);
end

figure; hold on
TT = bin_ctrs-1/num_bins/2;
for i = 1:count
    scatter(TT, SuccessRate_Binned(i, :), 6, 'm',  'filled')
end

X_plt = QT_Binned(:);
Y_plt = SuccessRate_Binned(:);
lm = fitlm(X_plt, Y_plt);
slope = lm.Coefficients.Estimate(2);
intercept = lm.Coefficients.Estimate(1);
fitted_values = predict(lm, X_plt);
plot(X_plt, fitted_values, 'k', 'LineWidth', 2);
R_squared = lm.Rsquared.Ordinary;
p_value = lm.Coefficients.pValue(2);
disp(['Slope: ', num2str(slope)]);
disp(['Intercept: ', num2str(intercept)]);
disp(['R-squared: ', num2str(R_squared)]);
disp(['p Value: ', num2str(p_value)]);


figure; hold on
ShadedPlot(QT_Binned(1, :), nanmean(SuccessRate_Binned, 1), 'k', 1, SEM(SuccessRate_Binned), [0.8 0.8 0.8])
plot(QT_Binned(1, :), nanmean(SuccessRate_Binned, 1), 'k', 'LineWidth', 1)
X_centroid = nanmean(QT_Binned, 1);
Y_centroid = nanmean(SuccessRate_Binned, 1);
Y_mean = nanmean(Y_centroid);
fitted_line = predict(lm, X_centroid');
SS_total = sum((Y_centroid'-Y_mean).^2);
SS_residual = sum((Y_centroid'-fitted_line).^2);
disp(['Fit mean R-squared: ', num2str(1-SS_residual/SS_total)]);


ylabel('Fraction of successful withheld trials')
xlabel('Mean $\beta$ weights', 'Interpreter', 'latex')
%%
figure; hold on
bar(1, nanmean(Slope_NE))
bar(2, nanmean(Slope_Ach))
bar(3, nanmean(Slope_sync))
errorbar([1 2 3], [nanmean(Slope_NE) nanmean(Slope_Ach) nanmean(Slope_sync)], ...
    [SEM(Slope_NE) SEM(Slope_Ach) SEM(Slope_sync)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Linear regression coefficient')
title('Success rate vs. mean $\beta$ weights', 'Interpreter', 'latex')




%%
%%
i = 26;
% Assuming 'A' is your 100-by-13 matrix
A = TRFcmp_NE.wtrans_coef_P{i};
B = TRFcmp_NE.wtrans_coef_R{i};
% Normalize the values within each row
validA = find(~isnan(A(:, 1)) & abs(mean(A(:, 4:end), 2))<=2);
% Normalize the values within each row
validB = find(~isnan(B(:, 1)) & abs(mean(B(:, 4:end), 2))<=2);

% validA = [];
% for i = 1:size(A, 1)
%     if isnan(A(i, 1))
%         continue
%     end
% %     if abs(mean(A(i, :))) > 2
% %         continue;
% %     end
%     if ~isempty(find(abs(A(i, :))>3))
%         continue;
%     end
%     validA(end+1) = i;
% end
% 
% validB = [];
% for i = 1:size(B, 1)
%     if isnan(B(i, 1))
%         continue
%     end
% %     if abs(mean(B(i, :))) > 2
% %         continue;
% %     end
%     if ~isempty(find(abs(B(i, :))>3))
%         continue;
%     end
%     validB(end+1) = i;
% end
%%
% figure; hold on
% scatter(zeros(1, size(A, 1)), nanmean(A(:, 4:end), 2))
% scatter(ones(1, size(B, 1)), nanmean(B(:, 4:end), 2))

%%
% A_normalized = UnitNormalization(A(validA, :)); % You can also use 'zscore' for z-score normalization, or your custom normalization function
A_normalized = A(validA, :);
figure
subplot(2, 2, 1)
% Create the heatmap
h = heatmap(A_normalized); % You can choose a different colormap if you prefer
colormap jet

% Define all x-axis labels
x_labels_all = {'-0.2', '', '0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'};

% Set the x-axis labels
h.XDisplayLabels = x_labels_all;

% Label the axes if needed
ylabel('Trials'); title('Failed trials')
h.YDisplayLabels = repmat({' '}, 1, length(h.YDisplayData));
% Add a colorbar to the heatmap to indicate the mapping between colors and values
h.GridVisible = 'off';
colorbar;
% h.Colorbar.Title.String = 'Weights a.u.';

% B_normalized = UnitNormalization(B(validB, :)); % You can also use 'zscore' for z-score normalization, or your custom normalization function
B_normalized = B(validB, :);
subplot(2, 2, 2)
% Create the heatmap
h = heatmap(B_normalized); % You can choose a different colormap if you prefer
colormap jet

% Define all x-axis labels
x_labels_all = {'-0.2', '', '0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'};

% Set the x-axis labels
h.XDisplayLabels = x_labels_all;

% Label the axes if needed
ylabel('Trials'); title('Successful trials')
h.YDisplayLabels = repmat({' '}, 1, length(h.YDisplayData));
% Add a colorbar to the heatmap to indicate the mapping between colors and values
h.GridVisible = 'off';
colorbar;
% h.Colorbar.Title.String = 'Weights a.u.';

subplot(2, 2, 3); hold on
ShadedPlot(TRFcmp_NE.model_fwd.t/1000, nanmean(A_normalized, 1),  'k', 1, SEM(A_normalized), [0.8 0.8 0.8])
plot(TRFcmp_NE.model_fwd.t/1000, nanmean(A_normalized, 1), '--k', 'LineWidth', 1.5)
xlabel('Time lag (s)'); ylabel('a.u.'); vline(0, ':k')

subplot(2, 2, 4); hold on
ShadedPlot(TRFcmp_NE.model_fwd.t/1000, nanmean(B_normalized, 1),  'k', 1, SEM(B_normalized), [0.8 0.8 0.8])
plot(TRFcmp_NE.model_fwd.t/1000, nanmean(B_normalized, 1), '-k', 'LineWidth', 1.5)
xlabel('Time lag (s)'); ylabel('a.u.'); vline(0, ':k')

    
%%
x_coord = 0;  % Adjust this value to your desired position