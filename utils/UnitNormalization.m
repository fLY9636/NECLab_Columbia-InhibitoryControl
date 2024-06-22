function Y = UnitNormalization(X)
% UnitNormalization normalizes values in the vector X into [0, 1]
if min(size(X))==1 
    Y = (X-min(X))/(max(X)-min(X));
else
    Y = [];
    for i = 1:size(X, 1)
        Y(end+1, :) = (X(i, :)-min(X(i, :)))/(max(X(i, :))-min(X(i, :)));
    end
end
end