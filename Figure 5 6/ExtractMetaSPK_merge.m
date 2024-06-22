%% Load behavior files and table
HOME = 'G:\My Drive\Response Inhibition Project\';
LAB1 = 'F:\Backup\BMEN_9100\'; LAB2 = 'M:\My Drive\Response Inhibition Project\'; LAB3 = 'G:\My Drive\Response Inhibition Project\';
ROOTDIR = LAB3;
Behavior_dir = [ROOTDIR '2021-22_Attention\NP 2023-12\IHB\Behavior\']; Behavior_files = natsortfiles(filename_scan(Behavior_dir));
% load([ROOTDIR '2021-22_Attention\NP 2023-12\MetaSPK_test_1'])
%% Initialization
INPUT = MetaSPK_test_1;
COL_ID = INPUT.ID;
COL_Date = INPUT.Date;
COL_SpkUnits = INPUT.Spike_unit;
COL_SpkTimes = INPUT.Spike_time;
COL_UnitMap = INPUT.Unit_map;
COL_UnitType = INPUT.Unit_type;
COL_UnitLevel = INPUT.Unit_level;
COL_Sync = INPUT.Pupil_in_imec;
COL_PTime = INPUT.Pupil_time;
%%
% Load NP data
folder_dir = 'D:\NP\';
animal_name = 'N4';
session_name = '2024_03_26_g0';
folder_name = [animal_name '_' session_name '\' animal_name '_' session_name '_imec0'];

spk_times = readNPY([folder_dir '\' folder_name '\' 'spike_times.npy']);
spk_clusters = readNPY([folder_dir '\' folder_name '\' 'spike_clusters.npy']);
cluster_group = readtable([folder_dir '\' folder_name '\' 'cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');
% cluster_group = cluster_group;

cluster_group.group{1}
%%
spk_temp = readNPY([folder_dir '\' folder_name '\' 'spike_templates.npy']);
temp = readNPY([folder_dir '\' folder_name '\' 'templates.npy']);
amp = readNPY([folder_dir '\' folder_name '\' 'amplitudes.npy']);

% Append new session to the table
load([folder_dir '\' folder_name '\' 'sync_channel.mat'])
ss = 71;
load(Behavior_files{ss}) % change
[onset1, score] = ExtractSYNC(vector, Data);

COL_ID{end+1, 1} = animal_name;
COL_Date{end+1, 1} = session_name;
COL_SpkUnits{end+1, 1} = spk_clusters;
COL_SpkTimes{end+1, 1} = spk_times;
COL_UnitMap{end+1, 1} = [cluster_group.cluster_id cluster_group.ch];

session_UnitType = zeros(height(cluster_group), 1);
for i = 1:height(cluster_group)
    if strcmp(cluster_group.KSLabel{i}, 'good')==1
        if strcmp(cluster_group.group{i}, 'mua')==1 || strcmp(cluster_group.group{i}, 'noise')==1
            session_UnitType(i) = 0;
        else
            session_UnitType(i) = 10;
        end
    else
        if strcmp(cluster_group.group{i}, 'good')==1
            session_UnitType(i) = 1;
        else
            session_UnitType(i) = 0;
        end
    end
end        
COL_UnitType{end+1, 1} = session_UnitType;

COL_UnitLevel{end+1, 1} = cluster_group.fr;
COL_Sync{end+1, 1} = onset1(1, :)';
COL_PTime{end+1, 1} = onset1(2, :)';

MetaSPK_test_1 = table(COL_ID, COL_Date, COL_SpkUnits, COL_SpkTimes, COL_UnitMap, COL_UnitType, COL_UnitLevel, COL_Sync, COL_PTime, ...
    'VariableNames', {'ID', 'Date', 'Spike_unit', 'Spike_time', 'Unit_map', 'Unit_type', 'Unit_level', 'Pupil_in_imec', 'Pupil_time'});

%% Save 
MetaSPK_test_1 = MetaSPK_test_1; save([ROOTDIR '2021-22_Attention\NP 2023-12\MetaSPK_test_1.mat'], 'MetaSPK_test_1', '-v7.3');


%% Create matrices for UnitMatch
unique_clusters = cluster_group.cluster_id(find(strcmp(cluster_group.KSLabel, 'good')==1));
numClusters = numel(unique_clusters);
[dummmy, tempSize2, tempSize3] = size(temp(1, :, :)); % Example, adjust as necessary
unit_sum_waveforms1 = zeros(numClusters, tempSize2, tempSize3);
unit_sum_waveforms2 = zeros(numClusters, tempSize2, tempSize3);
% Iterate over each unique cluster
clusterCounter = 1; % Index to keep track of the current cluster's results
for cluster_id = unique_clusters(:).'
    disp(cluster_id) % Display the current cluster_id
    % Find indices of spikes belonging to the current cluster
    cluster_indices = find(spk_clusters == cluster_id);
    % Extract the waveforms for spikes in the current cluster
    cluster_waveforms = temp(cluster_id + 1, :, :);
    % Calculate the midpoint of spk_clusters
    midpoint = round(length(spk_clusters) / 2);
    % Logical indices for the first and second halves
    first_half_indices = cluster_indices <= midpoint;
    second_half_indices = cluster_indices > midpoint;
    % Pre-allocate sum_waveforms for efficiency
    sum_waveforms1 = zeros([sum(first_half_indices), tempSize2, tempSize3]);
    sum_waveforms2 = zeros([sum(second_half_indices), tempSize2, tempSize3]);
    % Apply amp scaling for the first half
    if any(first_half_indices)
        sum_waveforms1 = repmat(cluster_waveforms, [sum(first_half_indices), 1, 1]) .* reshape(amp(cluster_indices(first_half_indices)), [sum(first_half_indices), 1, 1]);
    end
    % Apply amp scaling for the second half
    if any(second_half_indices)
        sum_waveforms2 = repmat(cluster_waveforms, [sum(second_half_indices), 1, 1]) .* reshape(amp(cluster_indices(second_half_indices)), [sum(second_half_indices), 1, 1]);
    end
    % Compute the mean across the first dimension
    unit_sum_waveforms1(clusterCounter, :, :) = nanmean(sum_waveforms1, 1);
    unit_sum_waveforms2(clusterCounter, :, :) = nanmean(sum_waveforms2, 1);
    clusterCounter = clusterCounter + 1; % Move to the next cluster's results
end
%%
numFiles = size(unit_sum_waveforms1, 1);
% mkdir([ROOTDIR '2021-22_Attention\NP 2023-12\UnitMatch feed\' animal_name '_' session_name '\RawWaveforms']);
mkdir([folder_dir '\' animal_name '_' session_name '\RawWaveforms']);

for i = 1:numFiles
    combinedSlice = cat(1, unit_sum_waveforms1(i, :, :), unit_sum_waveforms2(i, :, :));
    combinedSlice = permute(combinedSlice, [2, 3, 1]);
    size(combinedSlice)

    % Save each combined slice to a separate NumPy file
    filename = sprintf('unit_%d_%s.npy', i, session_name);
%     writeNPY(combinedSlice, [ROOTDIR '2021-22_Attention\NP 2023-12\UnitMatch feed\' animal_name '_' session_name '\RawWaveforms\' filename]);
    writeNPY(combinedSlice, [folder_dir '\' animal_name '_' session_name '\RawWaveforms\' filename]);
end
