function [h, p, t] = qw_statPairedTest(x,y)
% statistical test using either t-test or Mann-Whitney U test based on if x/y are from normal distribution
% t=1, t-test
% t=2, Mann-Whitney U test

if min(size(x))>1 && min(size(y))>1
    [H1,P1] = kstest((x-mean(x))/std(x));
    [H2,P2] = kstest((y-mean(y))/std(y));
    if H1 == 0 & H2 == 0
        [h p] = ttest(x, y);
        t = 1;
    else
        [p h] = signrank(x,y);
        t=2;
    end
else
    t = 1;
    [h p] = ttest(x, y);
end


