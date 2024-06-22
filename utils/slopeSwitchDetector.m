function numSwitches = slopeSwitchDetector(data, fs, windowSize)
    % Input:
    % data       - The data for which the slope is to be calculated
    % fs         - Sampling rate
    % windowSize - The size of the window over which linear regression is performed
    
    % Output:
    % numSwitches - Number of times the slope of the linear regression switches sign

    if length(data) < windowSize
        error('Data length must be greater than or equal to window size.');
    end

    numWindows = length(data) - windowSize + 1;
    slopes = zeros(1, numWindows);

    % Calculate the slopes for each window
    for i = 1:numWindows
        x = (0:windowSize-1)'/fs; % Time samples for current window
        y = data(i:i+windowSize-1);
        p = polyfit(x, y, 1); % First order polynomial fit (linear regression)
        slopes(i) = p(1); % Store the slope
    end

    % Check for sign changes in the slopes
    signs = sign(slopes);
    signDiffs = diff(signs);

    % Count the number of sign switches
    numSwitches = sum(abs(signDiffs) == 2);
end