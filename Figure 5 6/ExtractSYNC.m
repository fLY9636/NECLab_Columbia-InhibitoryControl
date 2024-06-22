function [onset1, score] = ExtractSYNC(vector, Data)
% ExtractSYNC has the inputs of vector (.mat file saving of imec-received
% TTLs from XPC) and Data (.mat file saving of XPC-sent TTLs to
% camera/imec). This function outputs the extracted rising-edges
% of received TTLs as indices in the imec timestamp.
edges = diff(vector);
onset1 = NaN(1, length(find(Data.Triggers.data(:, 2)==1)));                                                        % initialize the output     
score = 1; 
onset1(end+1, :) = Data.Triggers.data(find(Data.Triggers.data(:, 2)==1), end);                  % add cameratrig timestamp
onset = find(edges>=55 & edges<=65);
offset = find(edges>=-65 & edges<=-55);
disp(['Received: ' num2str(length(onset)) ' (RS)  ' num2str(length(offset)) ' (FL)  '])
disp(['Send: ' num2str(length(find(Data.Triggers.data(:, 2)==1)))])
if length(onset) < length(find(Data.Triggers.data(:, 2)))                                                       % when imec misses one or multiple triggers sent by XPC
    disp('receiver missing triggers')
    if isempty(find(diff(onset)>3600))                                                                                          % when the inter-trigger intervals are the same
        disp('missing the first/last trigger!')
        score = 0;
    else                                                                                                                                        % when missing triggers happened in the middle
        onset_interval = diff(onset);
        missing_after_idx = find(onset_interval>3600);                                                  
        fix_start_onset1 = 1;                                                                                                     
        fix_start_onset = 1;
        for i = 1:length(missing_after_idx)  % fix each missing trigger
            fix_end_onset = missing_after_idx(i);
            fix_end_onset1 = fix_start_onset1+fix_end_onset-fix_start_onset;
            onset1(1, fix_start_onset1:fix_end_onset1) = onset(fix_start_onset:fix_end_onset);
            fix_start_onset1 = fix_end_onset1+round(onset_interval(missing_after_idx(i))/3000);
            fix_start_onset = missing_after_idx(i)+1;
        end   
        onset1(1, fix_start_onset1:end) = onset(fix_start_onset:end);
        score = 1;
        disp('missing triggers labeled as NaN!')
        disp(['Min TTL itv: ' num2str(min(diff(onset1(1, :))))]) 
        disp(['Max TTL itv: ' num2str(max(diff(onset1(1, :))))]) 
    end
elseif length(onset) > length(find(Data.Triggers.data(:, 2)))                                                   % when imec receives more 'triggers' than sent by XPC
    disp('receiver has noise')
    error('please reset threshold!')
else                                                                                                                                          % when the triggers are equal on both sides
    disp('done')
    disp(['Min TTL itv: ' num2str(min(diff(onset)))])                                                                   % sanity check just in case
    disp(['Max TTL itv: ' num2str(max(diff(onset)))])                                                                   % sanity check just in case
    onset1(1, :) = onset;
    score = 1;
end
end

