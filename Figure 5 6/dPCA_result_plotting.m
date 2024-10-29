%% Plot average dPC traces across all sessions
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 207, 0);
loopIDX = 1:3;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
OFCIDX = 1:length(Behavior_files);
len = length(OFCIDX);

% IDX = [subjIDX{1}(10:end) subjIDX{2}(11:end) subjIDX{3}(9:end)]; % -- All no CNO session indices
IDX = [subjIDX{1}(end-4:end) subjIDX{2}(end-4:end) subjIDX{3}(end-4:end)]; % -- Saline session indices
% IDX = [subjIDX{1}(1:5) subjIDX{2}(1:5) subjIDX{3}(1:5)]; % -- Learning session indices
% IDX = [subjIDX{1}(1:end) subjIDX{2}(1:end) subjIDX{3}(1:end)]; % -- CNO session indices
% each cell represents a lump of one of three dPCs across sessions
LumpedSession_dpca_decision_Punish = cell(1, 3);
LumpedSession_dpca_decision_Success = cell(1, 3);
LumpedSession_dpca_independent_Punish = cell(1, 3);
LumpedSession_dpca_independent_Success = cell(1, 3);
for i = IDX(:).'
    for j = 1:3
        LumpedSession_dpca_decision_Punish{j}(end+1, :) = Lumped_dpca_decision{i, j}(1, :);
        LumpedSession_dpca_decision_Success{j}(end+1, :) = Lumped_dpca_decision{i, j}(2, :);
        LumpedSession_dpca_independent_Punish{j}(end+1, :) = Lumped_dpca_independent{i, j}(1, :);
        LumpedSession_dpca_independent_Success{j}(end+1, :) = Lumped_dpca_independent{i, j}(2, :);
    end
end
    
figure; hold on
for i = 1:3
    subplot(1, 3, i); hold on
    ShadedPlot(options.time, nanmean(LumpedSession_dpca_decision_Punish{i}, 1), 'k', 2, SEM(LumpedSession_dpca_decision_Punish{i}), [0.9 0.8 0.7])
    ShadedPlot(options.time, nanmean(LumpedSession_dpca_decision_Success{i}, 1), 'k', 2, SEM(LumpedSession_dpca_decision_Success{i}), [0.9 0.8 0.7])
    plot(options.time, nanmean(LumpedSession_dpca_decision_Punish{i}, 1), 'Color', 'm', 'LineWidth', 1.5)
    plot(options.time, nanmean(LumpedSession_dpca_decision_Success{i}, 1), 'Color', 'g', 'LineWidth', 1.5)
    title(['dPC ' num2str(i)])
    if i == 2
        xlabel('Time from behavioral outcome (s)')
    end
    if i == 1
        ylabel('Normalized firing rate (Hz)')
    end
end

figure; hold on
for i = 1:3
    subplot(1, 3, i); hold on
    ShadedPlot(options.time, nanmean(LumpedSession_dpca_independent_Punish{i}, 1), 'k', 2, SEM(LumpedSession_dpca_independent_Punish{i}), [0.9 0.8 0.7])
    ShadedPlot(options.time, nanmean(LumpedSession_dpca_independent_Success{i}, 1), 'k', 2, SEM(LumpedSession_dpca_independent_Success{i}), [0.9 0.8 0.7])
    plot(options.time, nanmean(LumpedSession_dpca_independent_Punish{i}, 1), 'Color', 'm', 'LineWidth', 1.5)
    plot(options.time, nanmean(LumpedSession_dpca_independent_Success{i}, 1), 'Color', 'g', 'LineWidth', 1.5)
    title(['dPC ' num2str(i)])
    if i == 2
        xlabel('Time from behavioral outcome (s)')
    end
    if i == 1
        ylabel('Normalized firing rate (Hz)')
    end
end

%% Plot dPC traces for each session
figure; hold on
for i = 1:length(IDX)
    subplot(5, 10, i); hold on
    title([num2str(round(table4Qi.SuccessRate(IDX(i)), 1)) ' %'])
    plot(options.time, LumpedSession_dpca_decision_Punish{2}(i, :), 'Color', 'm', 'LineWidth', 1)
    plot(options.time, LumpedSession_dpca_decision_Success{2}(i, :), 'Color', 'g', 'LineWidth', 1)
end

figure; hold on
for i = 1:length(IDX)
    subplot(5, 10, i); hold on
    title([num2str(round(table4Qi.SuccessRate(IDX(i)), 1)) ' %'])
    plot(options.time, LumpedSession_dpca_independent_Punish{1}(i, :), 'Color', 'm', 'LineWidth', 1)
    plot(options.time, LumpedSession_dpca_independent_Success{1}(i, :), 'Color', 'g', 'LineWidth', 1)
end

%% dPCs trajectory visualization (sessions)
% visualize in 2D space
figure; hold on
for i = 1:length(IDX)
%     plot(smoothdata(LumpedSession_dpca_decision_Punish{1}(i, :), 'movmean', 10), smoothdata(LumpedSession_dpca_decision_Punish{2}(i, :), 'movmean', 10), 'Color', [0.934 0.133 0.554 0.2], 'LineWidth', 2)
%     plot(smoothdata(LumpedSession_dpca_decision_Success{1}(i, :), 'movmean', 10), smoothdata(LumpedSession_dpca_decision_Success{2}(i, :), 'movmean', 10), 'Color', [0.015 0.861 0.089 0.2], 'LineWidth', 2)
%     scatter(LumpedSession_dpca_decision_Punish{1}(i, 1), LumpedSession_dpca_decision_Punish{2}(i, 1), 32, 'k', 'o', 'filled')
%     scatter(LumpedSession_dpca_decision_Punish{1}(i, end), LumpedSession_dpca_decision_Punish{2}(i, end), 32, 'k', 'x')
%     scatter(LumpedSession_dpca_decision_Success{1}(i, 1), LumpedSession_dpca_decision_Success{2}(i, 1), 32, 'k', 'o', 'filled')
%     scatter(LumpedSession_dpca_decision_Success{1}(i, end), LumpedSession_dpca_decision_Success{2}(i, end), 32, 'k', 'x')
    scatter(nanmean(LumpedSession_dpca_decision_Punish{1}(i, 1:80)), nanmean(LumpedSession_dpca_decision_Punish{2}(i, 1:80)), 32, [0.934 0.133 0.554], 'o', 'filled')
    scatter(nanmean(LumpedSession_dpca_decision_Success{1}(i, 1:80)), nanmean(LumpedSession_dpca_decision_Success{2}(i, 1:80)), 32, [0.015 0.861 0.089], 'o', 'filled')
end

xlabel('Demixed PC 1'); ylabel('Demixed PC 2');
%%
% visualize in 3D space
figure; hold on
for i = 1:length(IDX)
    % Plot for Punish in 3D
    plot3(smoothdata(LumpedSession_dpca_decision_Punish{1}(i, :), 'movmean', 5), ...
          smoothdata(LumpedSession_dpca_decision_Punish{2}(i, :), 'movmean', 5), ...
          smoothdata(LumpedSession_dpca_decision_Punish{3}(i, :), 'movmean', 5), ...
          'Color', [0.934 0.133 0.554], 'LineWidth', 0.5);
    
    % Plot for Success in 3D
    plot3(smoothdata(LumpedSession_dpca_decision_Success{1}(i, :), 'movmean', 5), ...
          smoothdata(LumpedSession_dpca_decision_Success{2}(i, :), 'movmean', 5), ...
          smoothdata(LumpedSession_dpca_decision_Success{3}(i, :), 'movmean', 5), ...
          'Color', [0.015 0.861 0.089], 'LineWidth', 0.5);
end
% for i = 1:length(IDX)
%     % Plot for Punish in 3D
%     plot3(LumpedSession_dpca_decision_Punish{1}(i, :), LumpedSession_dpca_decision_Punish{2}(i, :), LumpedSession_dpca_decision_Punish{3}(i, :), ...
%           'Color', [0.934 0.133 0.554], 'LineWidth', 0.5);
%     
%     % Plot for Success in 3D
%     plot3(LumpedSession_dpca_decision_Success{1}(i, :), LumpedSession_dpca_decision_Success{2}(i, :), LumpedSession_dpca_decision_Success{3}(i, :), ...
%           'Color', [0.015 0.861 0.089], 'LineWidth', 0.5);
% end
xlabel('dPC1'); ylabel('dPC2'); zlabel('dPC3');
grid on
view(3) % Set the default view to 3D
%% 2-D Euclidean distance
Dist_dpca = zeros(length(IDX), 1);
figure; hold on
for i = 1:length(IDX)
    subplot(5, 10, i); hold on
    title([num2str(round(table4Qi.NormSuccessRate(IDX(i)), 1)) ' %'])
    plot(LumpedSession_dpca_decision_Punish{1}(i, :), LumpedSession_dpca_decision_Punish{2}(i, :), 'Color', 'm', 'LineWidth', 1)
    scatter(LumpedSession_dpca_decision_Punish{1}(i, 1), LumpedSession_dpca_decision_Punish{2}(i, 1), 12, [0.85 0.6, 0.91], 'o', 'filled')
    scatter(LumpedSession_dpca_decision_Punish{1}(i, end), LumpedSession_dpca_decision_Punish{2}(i, end), 32, [0.85 0.6, 0.91], 'o', 'filled')
    plot(LumpedSession_dpca_decision_Success{1}(i, :), LumpedSession_dpca_decision_Success{2}(i, :), 'Color', 'g', 'LineWidth', 1)
    scatter(LumpedSession_dpca_decision_Success{1}(i, 1), LumpedSession_dpca_decision_Success{2}(i, 1), 12, [0.4 0.52, 1], 'o', 'filled')
    scatter(LumpedSession_dpca_decision_Success{1}(i, end), LumpedSession_dpca_decision_Success{2}(i, end), 32, [0.4 0.52, 1], 'o', 'filled')
    xlim([-10 10]); ylim([-10 10])
    
    Dist_dpca(i) = nanmean(sqrt((LumpedSession_dpca_decision_Punish{1}(i, 1:80)-LumpedSession_dpca_decision_Success{1}(i, 1:80)).^2+ ...
        (LumpedSession_dpca_decision_Punish{2}(i, 1:80)-LumpedSession_dpca_decision_Success{2}(i, 1:80)).^2));
end

%% Mahalanobis distance 2D
mahalanobisDistances = zeros(size(IDX));
TIDX = 1:80; % -5 tp -1s from OC
for i = 1:length(IDX)
    % Extract session data for Punish and Success
    punishData = [LumpedSession_dpca_decision_Punish{1}(i, TIDX); LumpedSession_dpca_decision_Punish{2}(i, TIDX)]';
    successData = [LumpedSession_dpca_decision_Success{1}(i, TIDX); LumpedSession_dpca_decision_Success{2}(i, TIDX)]'; % dPC1 and dPC2 during [-5 -1]s OC
    centroidPunish = mean(punishData);
    centroidSuccess = mean(successData); %centroids
    combinedData = [punishData; successData]; % covariance calculation
    covMatrix = cov(combinedData);
    difference = centroidPunish - centroidSuccess;
    mahalanobisDistances(i) = sqrt(difference / covMatrix * difference'); % Calculate Mahalanobis distance between the centroids
end
%% Mahalanobis distance 3D
mahalanobisDistances = zeros(size(IDX));
TIDX = 1:160; % -5 tp -1s from OC
for i = 1:length(IDX)
    % Extract session data for Punish and Success
%     punishData = [LumpedSession_dpca_decision_Punish{1}(i, TIDX); LumpedSession_dpca_decision_Punish{2}(i, TIDX); LumpedSession_dpca_decision_Punish{3}(i, TIDX)]';
%     successData = [LumpedSession_dpca_decision_Success{1}(i, TIDX); LumpedSession_dpca_decision_Success{2}(i, TIDX); LumpedSession_dpca_decision_Success{3}(i, TIDX)]'; % dPC1 and dPC2 during [-5 -1]s OC
    %************ movmean smoothed result ****************
    punishData = [smoothdata(LumpedSession_dpca_decision_Punish{1}(i, TIDX), 'movmean', 3); ...
        smoothdata(LumpedSession_dpca_decision_Punish{2}(i, TIDX), 'movmean', 3); ...
        smoothdata(LumpedSession_dpca_decision_Punish{3}(i, TIDX), 'movmean', 3)]';
    successData = [smoothdata(LumpedSession_dpca_decision_Success{1}(i, TIDX), 'movmean', 3); ...
        smoothdata(LumpedSession_dpca_decision_Success{2}(i, TIDX), 'movmean', 3); ...
        smoothdata(LumpedSession_dpca_decision_Success{3}(i, TIDX), 'movmean', 3)]'; % dPC1 and dPC2 during [-5 -1]s OC
    %********************************************************
    centroidPunish = mean(punishData);
    centroidSuccess = mean(successData); %centroids
    combinedData = [punishData; successData]; % covariance calculation
    covMatrix = cov(combinedData);
    difference = centroidPunish - centroidSuccess;
    mahalanobisDistances(i) = sqrt(difference / covMatrix * difference'); % Calculate Mahalanobis distance between the centroids
end

%%
IsolationIDX_Saline = mahalanobisDistances;
%%
IsolationIDX_CNO = mahalanobisDistances;
%% Plot linear regression fit between behavioral performance and dPC variables
figure; hold on
X = table4Qi.NormSuccessRate(IDX);
% Y = abs(nanmean(LumpedSession_dpca_decision_Punish{1}(:, 1:80)-LumpedSession_dpca_decision_Success{1}(:, 1:80), 2));
Y = IsolationIDX_CNO';
% Y = Dist_dpca;
valid_indices = isfinite(X) & isfinite(Y);
X_valid = X(valid_indices);
Y_valid = Y(valid_indices);
coefficients = polyfit(X_valid, Y_valid, 1);
slope = coefficients(1);
intercept = coefficients(2);
fitted_line = slope*X_valid+intercept;
scatter(X, Y)
scatter(X(1:5), Y(1:5), 32, 'b', 'o')
scatter(X(6:10), Y(6:10), 32, 'r', 'o')
scatter(X(11:15), Y(11:15), 32, 'g', 'o')
plot(X_valid, fitted_line, 'r', 'LineWidth', 2)
% Y_mean = mean(Y_valid);
% SS_total = sum((Y_valid-Y_mean).^2);
% SS_residual = sum((Y_valid-fitted_line).^2);
% R_squared = 1-SS_residual/SS_total;
% [p, ~, mu] = polyfit(X_valid, Y_valid, 1);
% disp(['Slope: ', num2str(slope)]);
% disp(['Intercept: ', num2str(intercept)]);
% disp(['R-squared: ', num2str(R_squared)]);
% disp(['p Value: ', num2str(p)]);
mdl = fitlm(X_valid, Y_valid)
xlabel('Normalized behavioral performance %')
ylabel('Cluster distance')

%% Plot average explained variance across conditions across sessions (generate Lumped variances from each condition first!)
figure; hold on
bar([1 2 3], [nanmean(Lumped_expVar_learning, 1) nanmean(Lumped_expVar_saline, 1) nanmean(Lumped_expVar_cno, 1)])
scatter(ones(length(Lumped_expVar_learning), 1), Lumped_expVar_learning)
scatter(2*ones(length(Lumped_expVar_saline), 1), Lumped_expVar_saline)
scatter(3*ones(length(Lumped_expVar_cno), 1), Lumped_expVar_cno)
errorbar([1 2 3], [nanmean(Lumped_expVar_learning, 1) nanmean(Lumped_expVar_saline, 1) nanmean(Lumped_expVar_cno, 1)], [SEM(Lumped_expVar_learning) SEM(Lumped_expVar_saline) SEM(Lumped_expVar_cno)], ...
    'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);


%% Barplot for comparing Saline vs CNO isolation index 
% comparison across sessions
figure; hold on
bar([1 2], [nanmean(IsolationIDX_Saline) nanmean(IsolationIDX_CNO)])
scatter(ones(size(IsolationIDX_Saline)), IsolationIDX_Saline)
scatter(2*ones(size(IsolationIDX_CNO)), IsolationIDX_CNO)
errorbar([1 2], [nanmean(IsolationIDX_Saline) nanmean(IsolationIDX_CNO)], [SEM(IsolationIDX_Saline') SEM(IsolationIDX_CNO')], ...
    'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Cluster distance')
% comparison across subjects
subj_IsolationIDX_Saline = [nanmean(IsolationIDX_Saline(1:5)) nanmean(IsolationIDX_Saline(6:10)) nanmean(IsolationIDX_Saline(11:15))];
subj_IsolationIDX_CNO = [nanmean(IsolationIDX_CNO(1:5)) nanmean(IsolationIDX_CNO(6:10)) nanmean(IsolationIDX_CNO(11:15))];
figure; hold on
bar([1 2], [nanmean(subj_IsolationIDX_Saline) nanmean(subj_IsolationIDX_CNO)])
for i = 1:3
    plot([1 2], [subj_IsolationIDX_Saline(i) subj_IsolationIDX_CNO(i)], '-ok')
end
errorbar([1 2], [nanmean(subj_IsolationIDX_Saline) nanmean(subj_IsolationIDX_CNO)], [SEM(subj_IsolationIDX_Saline') SEM(subj_IsolationIDX_CNO')], ...
    'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Cluster distance')