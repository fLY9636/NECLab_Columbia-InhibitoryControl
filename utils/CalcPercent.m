function Y = CalcPercent(X, baseline)
% CalcPercent caluclates the percentage change dynamics of X. The baseline
% is chosen as defined by range, which sets the start and end element and
% calculate the mean. Then Y = (X-mean(X(range)))/mean(X(range))*100
Y = (X-baseline)./baseline*100;
return