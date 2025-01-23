filePath = 'YOUR_DIRECTORY.csv'; % load data 
loadedTable = readtable(filePath);

% replace VAR,TABLE, FIELDNAME with your own
VAR1 = TABLENAME.FIELDNAME1; 
VAR2 = TABLENAME.FIELDNAME2; 
%% main - dots plot
% if length(VAR1)==length(VAR2) meaning they are paired sessions or animals
figure; hold on
bar([1 2], [nanmean(VAR1) nanmean(VAR2)])
for i = 1:length(VAR1) 
    plot([1 2], [VAR1(i) VAR2(i)], '-ok')
end
errorbar([1 2], [nanmean(VAR1) nanmean(VAR2)], [SEM(VAR1) SEM(VAR2)], ...
    'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);

% if length(VAR1)~=length(VAR2) meaning they are unpaired
% figure; hold on
% bar([1 2], [nanmean(VAR1) nanmean(VAR2)])
% scatter(ones(size(VAR1)), VAR1)
% scatter(2*ones(size(VAR2)), VAR2)
% errorbar([1 2], [nanmean(VAR1) nanmean(VAR2)], [SEM(VAR1) SEM(VAR2)], ...
%     'r', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
