function pValue = performAnova(dataMatrix)
% dataMatrix has the shape #ofSamples-by-#ofGroups
    % Check if the data matrix is empty or not
    if isempty(dataMatrix)
        error('Data matrix is empty.');
    end

    % Prepare the data and groups for ANOVA
    [numSamples, numGroups] = size(dataMatrix);
    groups = [];
    data = [];
    
    for i = 1:numGroups
        currentGroupData = dataMatrix(:, i);
        currentGroupData = currentGroupData(~isnan(currentGroupData)); % Remove NaN values
        groups = [groups; repmat(i, length(currentGroupData), 1)];
        data = [data; currentGroupData];
    end
    
    % Perform ANOVA
    [p, ~, stats] = anova1(data, groups, 'off'); % 'off' turns off the ANOVA table display

    % Display ANOVA statistics (optional)
    fprintf('ANOVA Statistics:\n');
    disp(stats);
    
    % Return the p-value
    pValue = p;
end