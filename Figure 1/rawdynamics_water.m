figure; hold on
subjMEAN_NE = [];
subjMEAN_Ach = [];
for i = 1:7
    subjMEAN_NE(end+1, :) = nanmean(cell2mat(Lumped_SpoutAdaption_NE(subjIDX{i})), 1);
    subjMEAN_Ach(end+1, :) = nanmean(cell2mat(Lumped_SpoutAdaption_Ach(subjIDX{i})), 1);
end
subplot(2, 2, [1 3]); hold on
ShadedPlot((WINst:WINed)/120, nanmean(subjMEAN_NE, 1), [0 0 1], 1, SEM(subjMEAN_NE), [0.73 0.83 0.96])
plot((WINst:WINed)/120, nanmean(subjMEAN_NE, 1), 'b', 'LineWidth', 1)
ShadedPlot((WINst:WINed)/120, nanmean(subjMEAN_Ach, 1), [1 0 0], 1, SEM(subjMEAN_Ach), [0.9 0.8 0.7])
plot((WINst:WINed)/120, nanmean(subjMEAN_Ach, 1), 'r', 'LineWidth', 1)
title('0.001Hz HPF')
xlabel('Time after Water Delivery (s)'); ylabel('Z-Score'); vline(0, '-k')

% Plot the peak value
subplot(2, 2, 2); hold on
bar(1, mean(Peak_NE(:, 2)))
bar(2, mean(Peak_Ach(:, 2)))
for i = 1:length(subjIDX)
    plot(1:2, [Peak_NE(i, 2) Peak_Ach(i, 2)], '-ok')
end
errorbar(1:2, [mean(Peak_NE(:, 2)) mean(Peak_Ach(:, 2))], [SEM(Peak_NE(:, 2)) SEM(Peak_Ach(:, 2))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Post-Water peak value - baseline (z-score)')

% Plot the peak latency
subplot(2, 2, 4); hold on
bar(1, mean(Peak_NE(:, 1)))
bar(2, mean(Peak_Ach(:, 1)))
for i = 1:length(subjIDX)
    plot(1:2, [Peak_NE(i, 1) Peak_Ach(i, 1)], '-ok')
end
errorbar(1:2, [mean(Peak_NE(:, 1)) mean(Peak_Ach(:, 1))], [SEM(Peak_NE(:, 1)) SEM(Peak_Ach(:, 1))], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Post-Water peak latency (s)')

figure;
% Plot subplot for NE heatmap
subplot(1, 2, 1);
imagesc(ds_NE_water);
colormap('hot');
colorbar('off');
caxis([min(ds_NE_water(:)), max(ds_NE_water(:))]); 
xlabel('Time after water delivery (s)');
xticks([1, 10, 20, 30, 40, 50, 60, 70]);
xticklabels({'-1', '0', '1', '2', '3', '4', '5', '6'});
ylabel('Trial');
yticks([1, size(ds_NE_water, 1)]);
vline(10, '-g')
% Plot subplot for ACh heatmap
subplot(1, 2, 2);
imagesc(ds_ACh_water);
colormap('hot');
caxis([min(ds_ACh_water(:)), max(ds_ACh_water(:))]); 
xlabel('Time after water delivery (s)');
xticks([1, 10, 20, 30, 40, 50, 60, 70]);
xticklabels({'-1', '0', '1', '2', '3', '4', '5', '6'});
yticks([]);
yticklabels([]);
vline(10, '-g')