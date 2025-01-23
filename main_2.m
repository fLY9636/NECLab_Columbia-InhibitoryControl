filePath = 'YOUR_DIRECTORY.csv'; % load data 
loadedTable = readtable(filePath);

% replace VAR,TABLE, FIELDNAME with your own
VAR1 = TABLENAME.FIELDNAME1; 
VAR2 = TABLENAME.FIELDNAME2; 
%% main - linear regression
X = VAR1;
Y = VAR2;

valid_indices = isfinite(X) & isfinite(Y);
X_valid = X(valid_indices);
Y_valid = Y(valid_indices);
coefficients = polyfit(X_valid, Y_valid, 1);
slope = coefficients(1);
intercept = coefficients(2);
fitted_line = slope*X_valid+intercept;

scatter(X, Y) % plot fitted line surrounded by point clouds
plot(X_valid, fitted_line, 'r', 'LineWidth', 2)

mdl = fitlm(X_valid, Y_valid) % will generate stats