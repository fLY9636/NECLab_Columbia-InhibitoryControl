function PlotSpikeTrain(A, fs)
% h=PlotSpikeTrain(A) will plot the spike trains (COLUMN AXIS) across all trials (ROW AXIS)
% A has to be a matrix
yHeight = 0.8; % Define the height of each vertical line

% Convert the 0s and 1s in A to the desired y-coordinates
y = repmat((1:size(A, 1))', 1, size(A, 2)).*A*yHeight;

% Plot the y-coordinates using imagesc with a colormap that maps 1s to black and 0s to white
imagesc(y, [0 yHeight]);
colormap([1 1 1; 0 0 0]); % Use a colormap that maps 1s to black and 0s to white

time = (0:size(A, 2))/fs-0.1*fs;
xlabel('Time after Inhibition Tone Onset (s)')
ylabel('Trial')
xticks(linspace(1, size(A, 2), 5));
xticklabels(sprintf('%d\n', floor(time(xticks))));
end