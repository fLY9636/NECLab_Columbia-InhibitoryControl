function fig = ShowColorCode(subj_COLOR)
    % Calculate patch width and height
    num_colors = length(subj_COLOR);
    patch_width = 20;
    patch_height = 20;

    % Create figure
    fig = figure;
    hold on;

    % Plot color patches and add text
    for i = 1:num_colors
        % Calculate patch position
        x = (i - 1) * patch_width;
        y = 0;

        % Draw patch
        patch([x x+patch_width x+patch_width x], [y y y+patch_height y+patch_height], subj_COLOR{i}, 'EdgeColor', 'none');

        % Calculate text position (center of patch)
        text_x = x + patch_width / 2;
        text_y = y + patch_height / 2;

        % Add text
        text(text_x, text_y, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'w');
    end

    % Set axis properties
    xlim([0, num_colors * patch_width]);
    ylim([0, patch_height]);
    axis off;

    % Add title
    title('Color Band');

    % Hold off
    hold off;
end