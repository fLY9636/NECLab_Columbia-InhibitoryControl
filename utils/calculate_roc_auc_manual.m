function AUC = calculate_roc_auc_manual(vector1, vector2)
    % Combine the two vectors to create labels (0 for vector1, 1 for vector2)
    labels = [zeros(size(vector1)); ones(size(vector2))];
    
    % Combine the two vectors to create the corresponding scores
    scores = [vector1; vector2];
    
    % Sort the scores in descending order and their corresponding labels
    [sorted_scores, idx] = sort(scores, 'descend');
    sorted_labels = labels(idx);
    
    % Calculate the number of positive and negative instances
    num_positives = sum(sorted_labels == 1);
    num_negatives = sum(sorted_labels == 0);
    
    % Initialize variables to keep track of True Positive (TP) and False Positive (FP) counts
    TP_count = 0;
    FP_count = 0;
    
    % Initialize variables to keep track of the previous TPR and FPR
    prev_TPR = 0;
    prev_FPR = 0;
    
    % Initialize the area under the ROC curve (AUC) as zero
    AUC = 0;
    
    % Calculate the area under the ROC curve using the trapezoidal rule
    for i = 1:length(sorted_scores)
        if sorted_labels(i) == 1
            TP_count = TP_count + 1;
        else
            FP_count = FP_count + 1;
        end
        
        TPR = TP_count / num_positives;
        FPR = FP_count / num_negatives;
        
        AUC = AUC + (TPR + prev_TPR) * (FPR - prev_FPR) / 2;
        
        prev_TPR = TPR;
        prev_FPR = FPR;
    end
end
