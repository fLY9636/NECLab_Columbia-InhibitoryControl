function [MeanINBin, SEMINBin, bin_ctrs, bin_IDXs] = BinnedQTMean1d(X, Y, num_bins)
% BinnedMean function divides X into num_bins, extracts the elements which
% falls into each bin (Quantiles), and calculates the mean of corresponding elements in
% Y
MEAN = []; sem = []; bin_IDXs = {};
bin_ctrs = 0:1/num_bins:1; 
QT = quantile(X, bin_ctrs);
bin_ctrs = bin_ctrs(1:end-1)+1/num_bins;
for i = 1:num_bins
    IDXbin = find(X>=QT(i) & X<=QT(i+1)); bin_IDXs{end+1, 1} = IDXbin;
    MEAN(end+1, 1) = mean(Y(IDXbin));
    sem(end+1, 1) = SEM(Y(IDXbin));
end
MeanINBin = MEAN;
SEMINBin = sem;
end