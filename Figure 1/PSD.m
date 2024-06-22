figure; hold on
subjC_set = [];
cnt = 0;
% plot averaged coherence across subjects and individual coherence 
for subj = 1:length(subjIDX)
    sessionmeanC_set = [];
    cnt = cnt+1;
    for i = subjIDX{subj}(:).'
        cnt = cnt+1;
        sessionmeanC_set(end+1, :) = nanmean(reshape(C_Lumped(:, :, i), size(C_Lumped, 1), size(C_Lumped, 2)), 1);
    end
    subjC_set(end+1, :) = nanmean(sessionmeanC_set, 1);
%     plot(f1, nanmean(sessionmeanC_set, 1), 'k', 'LineWidth', 0.5)
end
idx = find(f1>=0.05);
ShadedPlot(f1(idx), nanmean(subjC_set(:, idx), 1), [0 0 1], 2, SEM(subjC_set(:, idx)), [0.8 0.9 0.9])
plot(f1, nanmean(subjC_set, 1), 'b', 'LineWidth', 1)
xlim([0.05 3.5])

pltIDX = find(f1>=0.05);
h = patch([f1(pltIDX) fliplr(f1(pltIDX))], [QL fliplr(QU)], [0, .6, .77], 'FaceAlpha',0.5, 'EdgeColor','none')