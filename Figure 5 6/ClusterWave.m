%% Plot example session units' average waveforms with SEM
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 200, 0);
loopIDX = 1:3;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
OFCIDX = 1:length(Behavior_files);
len = length(OFCIDX);

t_waveform = -41:40;
METAMATRIX = MetaSPK_test_1; 
figure; hold on
for i = 1
    load(Behavior_files{i})
    load([ROOTDIR '2021-22_Attention\NP 2023-12\waveforms\longer IHB\' 'Session_unit_waveform_' Data.AnimalID '_' num2str(Data.date(1)) '-' num2str(Data.date(2)) '-' num2str(Data.date(3)) '-' num2str(Data.date(4)) ...
            '.mat'])
    goodunit_idx = find(METAMATRIX.Unit_type{i}~=0 & METAMATRIX.Unit_level{i}>1);
    channel_ids = METAMATRIX.Unit_map{i}(goodunit_idx, 2); 
    subplot(1, 1, i); hold on
    for j = 1:size(waveFormsMeanSEM.MEAN, 1)
        ShadedPlot(t_waveform, squeeze(waveFormsMeanSEM.MEAN(j, channel_ids(j)+1, :))', 'k', 1, squeeze(waveFormsMeanSEM.SEM(j, channel_ids(j)+1, :))', [0.8 0.8 0.8])
        plot(t_waveform, squeeze(waveFormsMeanSEM.MEAN(j, channel_ids(j)+1, :)))
    end
end
%% Calculate average waveforms and extract h5 file for UMAP analysis

%% Calculate and collect average waveforms from each session (once saved average waveforms, NO NEED TO RUN THIS)
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 200, 0);
loopIDX = 1:3;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
OFCIDX = 1:length(Behavior_files);
len = length(OFCIDX);

IDX = [subjIDX{1}(end-4:end) subjIDX{2}(end-4:end) subjIDX{3}(end-4:end)]; 
data_dir = [ROOTDIR '2021-22_Attention\NP 2023-12\waveforms\longer IHB']; data_files = natsortfiles(filename_scan(data_dir)); % chosse Saline or DCZ dataset
t_waveform = -41:120;
METAMATRIX = MetaSPK_test_1;  % choose Saline or DCZ dataset
Lumped_average_waveform = cell(length(IDX), 1);
cnt = 0;
for i = IDX(:).'
    i
    cnt = cnt+1;
    load(data_files{i})
    goodunit_idx = find(METAMATRIX.Unit_type{i}~=0 & METAMATRIX.Unit_level{i}>1);
    channel_ids = METAMATRIX.Unit_map{i}(goodunit_idx, 2);
    for j = 1:size(waveFormsMeanSEM.MEAN, 1)
        Lumped_average_waveform{cnt}(end+1, :) = waveFormsMeanSEM.MEAN(j, channel_ids(j)+1, :); % for fast good units
%         Lumped_average_waveform{cnt}(end+1, :) = waveFormsMeanSEM.MEAN(j, 1, :); % for all good units
    end
end

%% Feature extraction for averaged waveform of each unit (AP width [half min],  Trough2Peak, PrePost Peak ratio)
Lumped_waveform_features = cell(length(Lumped_average_waveform), 1); 
for i = 1:length(Lumped_average_waveform)
    waveforms = Lumped_average_waveform{i}; % Get the waveforms for the current cell
    features_matrix = zeros(size(waveforms, 1), 3); % Preallocate matrix for features of all neurons
    for j = 1:size(waveforms, 1)
        waveform = zscore(waveforms(j, :));
        features_matrix(j, :) = ExtractWaveformFeatures_corr(waveform, 20);
    end
    Lumped_waveform_features{i} = features_matrix; % Assign features matrix to the corresponding cell
end

% Lumped_waveform_features_NORM = Lumped_waveform_features;
% cnt = 0;
% for i = 1:length(subjIDX)
%     NormPool = cell2mat(Lumped_waveform_features((i-1)*5+1:i*5)); % gather session units' features from the same subject
%     for j = 1:5
%         cnt = cnt+1;
%         Lumped_waveform_features_NORM{cnt}(:, 1) = (Lumped_waveform_features_NORM{cnt}(:, 1)-nanmean(NormPool(:, 1)))/nanstd(NormPool(:, 1));
%         Lumped_waveform_features_NORM{cnt}(:, 2) = (Lumped_waveform_features_NORM{cnt}(:, 2)-nanmean(NormPool(:, 2)))/nanstd(NormPool(:, 2));
%         Lumped_waveform_features_NORM{cnt}(:, 3) = (Lumped_waveform_features_NORM{cnt}(:, 3)-nanmean(NormPool(:, 3)))/nanstd(NormPool(:, 3));
%     end
% end

figure; hold on
MixedSubjFeatures = cell2mat(Lumped_waveform_features);
% MixedSubjFeatures = cell2mat(Lumped_waveform_features);
scatter(MixedSubjFeatures(:, 2), MixedSubjFeatures(:, 1))
% scatter(MixedSubjFeatures(:, 2), MixedSubjFeatures(:, 1))

% figure; hold on
% scatter3(MixedSubjFeatures(:, 1), MixedSubjFeatures(:, 2), MixedSubjFeatures(:, 3))

%% Feature extraction for averaged waveform of each unit (AP width [half max],  Trough2Peak)
Lumped_waveform_features = cell(length(Lumped_average_waveform), 1); 
for i = 1:length(Lumped_average_waveform)
    waveforms = Lumped_average_waveform{i}; % Get the waveforms for the current cell
    features_matrix = zeros(size(waveforms, 1), 2); % Preallocate matrix for features of all neurons
    for j = 1:size(waveforms, 1)
        waveform = zscore(waveforms(j, :));
        features_matrix(j, :) = ExtractWaveformFeatures_3_corr(waveform, 20); % use 2 features or 3 features or consider performing interpolation first 
    end
    Lumped_waveform_features{i} = features_matrix; % Assign features matrix to the corresponding cell
end

Lumped_waveform_features_NORM = Lumped_waveform_features;
cnt = 0;
for i = 1:length(subjIDX)
    NormPool = cell2mat(Lumped_waveform_features((i-1)*5+1:i*5)); % gather session units' features from the same subject
    for j = 1:5
        cnt = cnt+1;
        Lumped_waveform_features_NORM{cnt}(:, 1) = (Lumped_waveform_features_NORM{cnt}(:, 1)-nanmean(NormPool(:, 1)))/nanstd(NormPool(:, 1));
        Lumped_waveform_features_NORM{cnt}(:, 2) = (Lumped_waveform_features_NORM{cnt}(:, 2)-nanmean(NormPool(:, 2)))/nanstd(NormPool(:, 2));
    end
end

figure; hold on
% MixedSubjFeaturesNORM = cell2mat(Lumped_waveform_features_NORM);
MixedSubjFeatures = cell2mat(Lumped_waveform_features);
% scatter(MixedSubjFeaturesNORM(:, 2), MixedSubjFeaturesNORM(:, 1))
scatter(MixedSubjFeatures(:, 2), MixedSubjFeatures(:, 1))

% figure; hold on
% scatter3(MixedSubjFeaturesNORM(:, 1), MixedSubjFeaturesNORM(:, 2), MixedSubjFeaturesNORM(:, 3))


%% GMM clustering 2D
 
figure
subplot(1, 3, 1) 
k = 2; % Define the number of clusters
options = statset('MaxIter',1000); % Increase max iterations if needed

% Fit GMM to the data
gmm = fitgmdist(MixedSubjFeatures, k, 'Options', options, 'RegularizationValue', 0.01);
cluster_labels = cluster(gmm, MixedSubjFeatures);
scatter(MixedSubjFeatures(:,2)/30, MixedSubjFeatures(:,1)/30, 16, cluster_labels);
ylabel('AP Width (ms)');
xlabel('Trough to Peak Duration (ms)');
% title('3D Plot of Waveform Features Clustering');
colormap(jet(k)); % Color map with 'k' distinct colors
colorbar;

% Calculate the frequency of each cluster label
clusterCounts = histcounts(cluster_labels, 1:k+1); % Using k clusters

% Create a pie chart of the cluster distributions
subplot(1, 3, 2)
hPie = pie(clusterCounts);
% title('Cluster Distribution');

% Get the jet colormap used in the GMM scatter plot
colormap(jet(k));
cmap = colormap; % Store the colormap

% % Apply the colormap to the pie chart
% for i = 1:length(hPie)
%     if mod(i, 2) == 1 % Odd entries in hPie are the pie slices
%         set(hPie(i), 'FaceColor', cmap(i-mod(i-1, 2)/2, :));
%     end
% end

%% GMM clustering 3D
figure;

% Subplot for 3D scatter plot of GMM clustering results
subplot(1, 3, 1)
k = 2; % Define the number of clusters
options = statset('MaxIter', 1000); % Increase max iterations if needed

% Fit GMM to the data
gmm = fitgmdist(MixedSubjFeatures, k, 'Options', options, 'RegularizationValue', 0.01);
cluster_labels = cluster(gmm, MixedSubjFeatures);


% 3D scatter plot for 2D data with an artificial third dimension
scatter3(MixedSubjFeatures(:, 3)/30, MixedSubjFeatures(:,2)/30, MixedSubjFeatures(:,1)/30, 12, cluster_labels, 'filled');
zlabel('AP Width (ms)');
ylabel('Trough to Peak Duration (ms)');
xlabel('Peak Ratio');
set(gca, 'XDir', 'reverse');
% title('3D Plot of 2D Waveform Features Clustering');
colormap(jet(k)); % Color map with 'k' distinct colors
colorbar;

% figure; hold on
% Calculate the frequency of each cluster label
clusterCounts = histcounts(cluster_labels, 1:k+1); % Using k clusters

% Subplot for pie chart of cluster distributions
subplot(1, 3, 2)
hPie = pie(clusterCounts);
title('Cluster Distribution');

% Get the jet colormap used in the GMM scatter plot
colormap(jet(k));
cmap = colormap; % Store the colormap

%% Plot average waveforms across clusters
% figure; hold on
numClusters = max(cluster_labels);  % Number of clusters identified
timePoints = (-41:120)/30000*1000;  % Time points for the waveform data
clusteredWaveforms = cell(numClusters, 1);
all_waveforms = cell2mat(Lumped_average_waveform);

% Separate waveforms by clusters
for k = 1:numClusters
    clusteredWaveforms{k} = all_waveforms(cluster_labels == k, :);
end

% Initialize figure for plotting
subplot(1, 3, 3) % following GMM clustering
hold on;

colors = jet(k);  % Get distinct colors for each cluster

% Plot each cluster's waveforms mean and SEM
for k = 1:numClusters
    % Calculate mean and SEM of the waveforms in the cluster
    meanWaveform = mean(clusteredWaveforms{k}, 1);
    semWaveform = std(clusteredWaveforms{k}, 0, 1) / sqrt(size(clusteredWaveforms{k}, 1));
    
    % Define the time and waveform plots for mean ± SEM
    upperSEM = meanWaveform + semWaveform;
    lowerSEM = meanWaveform - semWaveform;

    % Plot with shaded SEM area
    fill([timePoints, fliplr(timePoints)], [upperSEM, fliplr(lowerSEM)], ...
        colors(k,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(timePoints, meanWaveform, 'Color', colors(k,:), 'LineWidth', 2);
end

xlabel('Time (ms)');
ylabel('Amplitude');
% title('Mean ± SEM Waveforms by Cluster');
legend(arrayfun(@(k) sprintf('Cluster %d', k), 1:numClusters, 'UniformOutput', false));
% grid on;
hold off;

figure
for k = 1:numClusters
    subplot(1, 2, k); hold on
    for j = 1:size(clusteredWaveforms{k}, 1)
        plot(timePoints, clusteredWaveforms{k}(j, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    end
    plot(timePoints, mean(clusteredWaveforms{k}, 1), 'Color', colors(k,:), 'LineWidth', 2);
end

%% Look at BIC
maxClusters = 9; % Maximum number of clusters to test
BIC = zeros(maxClusters, 1); % Preallocate array for BIC values

options = statset('MaxIter', 1000); % Setting options for GMM

% Loop over the number of clusters k
for k = 1:maxClusters
    % Fit Gaussian Mixture Model
    gmm = fitgmdist(MixedSubjFeatures, k, 'Options', options, 'RegularizationValue', 0.01);
    % Store the BIC value
    BIC(k) = gmm.BIC;
end

% Plot BIC values
figure;
plot(1:maxClusters, BIC, '-o');
xlabel('Number of Clusters');
ylabel('Bayesian Information Criterion');
% title('BIC vs. Number of Clusters');
% grid on;










    











