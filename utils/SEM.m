function y = SEM(x)
% SEM calculates the standard error mean of x
% if x is a matrix, then it calculates the SEM of each column, as each row
% represents a sample
sample_size = size(x, 1);
STD = nanstd(x, 0, 1);
y = STD/sqrt(sample_size);
end