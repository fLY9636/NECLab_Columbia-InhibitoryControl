function h = ShadedPlot(x, MEAN, Color1, width, STD, Color2)
% this will generate a plot with shaded area indicating the STD around
% MEAN, using the built-in Color codes
% NOTE: x and MEAN must be a row vector!
% h = plot(x, MEAN, 'Color', Color1, 'LineWidth', width)
% hold all;
fill([x fliplr(x)], [(MEAN-STD) fliplr(MEAN+STD)], Color2, 'linestyle', 'none', 'FaceAlpha', 0.4)
