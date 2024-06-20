function features = ExtractWaveformFeatures_1(waveform, interpolationFactor)
% 1-> AP Peak Width  2-> Trough to Peak Duration
    % Initialize features vector
    features = zeros(1, 2);
    
    % Time vector for the original waveform
    originalTime = 1:length(waveform);
    
    % New time vector for the interpolated waveform
    interpolatedTime = linspace(1, length(waveform), length(waveform) * interpolationFactor);
    
    % Interpolate waveform
    interpolatedWaveform = interp1(originalTime, waveform, interpolatedTime, 'spline');
    
    % 1. Calculate the AP Peak Width (Full-Width Half Maximum of the peak)
    peakIndex = find(interpolatedWaveform == max(interpolatedWaveform), 1); % Find the maximum point (peak)
    halfMax = max(interpolatedWaveform) / 2;
    leftIndex = find(interpolatedWaveform(1:peakIndex) >= halfMax, 1, 'first');
    rightIndex = peakIndex + find(interpolatedWaveform(peakIndex:end) >= halfMax, 1, 'last') - 1;
    if isempty(leftIndex) || isempty(rightIndex)
        features(1) = NaN;
    else
        features(1) = (rightIndex - leftIndex) / interpolationFactor;
    end

    % 2. Calculate Trough to Peak Duration
    troughIndex = find(interpolatedWaveform == min(interpolatedWaveform), 1); % Find the minimum point (trough)
    if peakIndex > troughIndex % Ensure the peak follows the trough
        features(2) = (peakIndex - troughIndex) / interpolationFactor;
    else
        features(2) = NaN; % If peak does not follow trough, set to NaN
    end
end