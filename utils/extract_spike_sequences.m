function outputMatrix = extract_spike_sequences(X, Xs, T, fs)
    % X: Time series column vector
    % Xs: Spike event time vector
    % T: Window length (sec)
    % fs: Sampling rate
    
    num_spikes = length(Xs);
    outputMatrix = zeros(num_spikes, T*fs+1);
    validRows = true(num_spikes, 1); % Track valid rows
    
    for i = 1:num_spikes
        spike_time = Xs(i);
        start_idx = max(1, spike_time-round((T/2)*fs));
        end_idx = min(length(X), spike_time+round((T/2)*fs));
        
        % Check if there are enough samples to form the sequence
        if end_idx-start_idx+1>=T*fs
            % Extract sequence and store it in the output matrix
            sequence = X(start_idx:end_idx)';
            outputMatrix(i, :) = sequence;
        else
            % Handle the case where there are not enough samples
            fprintf('Warning: Insufficient samples for spike event %d\n', i);
            validRows(i) = false; % Mark the row as invalid
        end
    end
    
    % Remove rows with insufficient samples
    outputMatrix = outputMatrix(validRows, :);
end
