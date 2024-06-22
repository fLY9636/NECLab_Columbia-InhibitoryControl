function [FR0, SDF0] = sdfEstimate(spiketrains, step_size, window_size0, n, stdev, fs)
% sdfEstimate calculates the spike density function SDF0 and firing rate FR0 of the N-by-T spiketrains using
% moving window of window_size0
% FR1 and SDF1 are the Gaussian-smoothed results using kernel of
% window_size1 and stdev
spiketrains1 = [];
for i = 1:size(spiketrains, 1)
    if isempty(find(isnan(spiketrains(i, :))))
        spiketrains1(end+1, :) = spiketrains(i, :);
    end
end
[N, T] = size(spiketrains1);
spikecnt_matrix = [];
TT = window_size0*fs/2+1:step_size*fs:T-window_size0*fs/2;
for i = window_size0*fs/2+1:step_size*fs:T-window_size0*fs/2 % e.g., moving window is [-50, 50] in the case of window_size has 100 samples
    spikecnt_matrix(:, end+1) = sum(spiketrains1(:, i-window_size0*fs/2:i+window_size0*fs/2), 2, 'omitnan');
end
FR0 = zeros(1, size(spikecnt_matrix, 2)); SDF0 = zeros(1, size(spikecnt_matrix, 2));
for i = 1:size(spikecnt_matrix, 2)
    FR0(i) = sum(spikecnt_matrix(:, i), 1, 'omitnan')/window_size0/N; % instantaneous firing rate
    SDF0(i) = length(find(spikecnt_matrix(:, i)>0))/N; % instantaneous spike density function
end
% FR1 = smoothts(FR0, 'g', n, stdev); % Gaussian-smoothed results
% SDF1 = smoothts(SDF0, 'g', n, stdev);
return
