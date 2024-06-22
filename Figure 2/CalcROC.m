function [points, AUC] = CalcROC(distr1, distr2, thresholds)
% CalcROC calculates the Hit-FA rate pairs (as in points) and the AUC,
% given two distributions distr1 and distr2 within its range, and a set of detection
% thresholds. distr1 and distr2 should be the distribution of true and
% false signal respectively
TPR_set = [];
FPR_set = [];
for i = 1:length(thresholds)
    threshold = thresholds(i);
    TP = length(find(distr1>=threshold));
    FP = length(find(distr2>=threshold));
    TPR_set(end+1) = TP/length(distr1);
    FPR_set(end+1) = FP/length(distr2);
end
points = [FPR_set; TPR_set];
AUC = trapz(sort(FPR_set, 'ascend'), sort(TPR_set, 'ascend'));
end