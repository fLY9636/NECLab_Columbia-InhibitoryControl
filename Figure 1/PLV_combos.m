% PLV ~->pupil
figure; hold on                                                                                                            
ShadedPlot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_NElockpupil, 1), [0 0 1], 1, SEM(Lumped_NElockpupil), [0.73 0.83 0.96])
plot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_NElockpupil, 1), 'Color', [0 0 1], 'LineWidth', 1)
ShadedPlot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_Achlockpupil, 1), [1 0 0], 1, SEM(Lumped_Achlockpupil), [0.9 0.8 0.7])
plot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_Achlockpupil, 1), 'Color', [1 0 0], 'LineWidth', 1)
plot(rad2deg(bin_ctrs(1:num_bins/2)*6.2832-3.1416), 0.8*cos(bin_ctrs(1:num_bins/2)*6.2832-3.1416), 'Color', [1.00 0.60 0.78], 'LineStyle', ':', 'LineWidth', 2)
plot(rad2deg(bin_ctrs(num_bins/2+1:num_bins)*6.2832-3.1416), 0.8*cos(bin_ctrs(num_bins/2+1:num_bins)*6.2832-3.1416), 'Color', [0.68 0.92 1.00], 'LineStyle', ':', 'LineWidth', 2)
vline(0, '--k')

% PLV ~->NE
figure; hold on                                                                                                            
ShadedPlot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_AchlockNE, 1), [1 0 0], 1, SEM(Lumped_AchlockNE), [0.9 0.8 0.7])
plot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_AchlockNE, 1), 'Color', [1 0 0], 'LineWidth', 1)
ShadedPlot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_pupillockNE, 1), [0 0 0], 1, SEM(Lumped_pupillockNE), [0.8 0.8 0.8])
plot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_pupillockNE, 1), 'Color', [0 0 0], 'LineWidth', 1)
plot(rad2deg(bin_ctrs(1:num_bins/2)*6.2832-3.1416), 0.8*cos(bin_ctrs(1:num_bins/2)*6.2832-3.1416), 'Color', [1.00 0.60 0.78], 'LineStyle', ':', 'LineWidth', 2)
plot(rad2deg(bin_ctrs(num_bins/2+1:num_bins)*6.2832-3.1416), 0.8*cos(bin_ctrs(num_bins/2+1:num_bins)*6.2832-3.1416), 'Color', [0.68 0.92 1.00], 'LineStyle', ':', 'LineWidth', 2)
vline(0, '--k')

% PLV ~->ACh
figure; hold on                                                                                                            
ShadedPlot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_NElockAch, 1), [0 0 1], 1, SEM(Lumped_NElockAch), [0.73 0.83 0.96])
plot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_NElockAch, 1), 'Color', [0 0 1], 'LineWidth', 1)
ShadedPlot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_pupillockAch, 1), [0 0 0], 1, SEM(Lumped_pupillockAch), [0.8 0.8 0.8])
plot(rad2deg(bin_ctrs*6.2832-3.1416), nanmean(Lumped_pupillockAch, 1), 'Color', [0 0 0], 'LineWidth', 1)
plot(rad2deg(bin_ctrs(1:num_bins/2)*6.2832-3.1416), 0.8*cos(bin_ctrs(1:num_bins/2)*6.2832-3.1416), 'Color', [1.00 0.60 0.78], 'LineStyle', ':', 'LineWidth', 2)
plot(rad2deg(bin_ctrs(num_bins/2+1:num_bins)*6.2832-3.1416), 0.8*cos(bin_ctrs(num_bins/2+1:num_bins)*6.2832-3.1416), 'Color', [0.68 0.92 1.00], 'LineStyle', ':', 'LineWidth', 2)
vline(0, '--k')



