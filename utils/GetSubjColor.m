% Define the number of colors to generate
num_colors = 28;

% Generate random RGB values
rgb_values = rand(num_colors,3);

% Scale the RGB values to the range [0, 255]
rgb_values = rgb_values;

% Display the RGB values
disp(rgb_values);

% Convert RGB values to 1-by-3 vectors
color_vectors = num2cell(rgb_values, 2);

% Display the color vectors
disp(color_vectors);

subj_COLOR = color_vectors;