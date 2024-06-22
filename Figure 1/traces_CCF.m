



figure
% Plot sensors fluctuation
subplot(2, 1, 1); hold on
h1 = plot(time_axis, zraw_NE_fluc, '-b', 'LineWidth', 1.5)
h2 = plot(time_axis, zraw_ACh_fluc, '-r', 'LineWidth', 1.5)
legend([h1 h2], {'NE', 'ACh'})
xlabel('Time (s)')
ylabel('Z-score')
% Plot pupil fluctuation
subplot(2, 1, 2); 
plot(time_axis, zraw_pupil_fluc, '-k', 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Z-score')