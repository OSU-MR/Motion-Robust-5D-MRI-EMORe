function entries_below_const = plot_out(max_likelihood,const,fig_num,k,str,folder_out_path)

% Calculate the number of entries below 'const'

% Plot the histogram
figure(fig_num);
histogram(max_likelihood(:)); % 'pdf' for probability density function
hold on; % Allow overlaying plots
if ~isempty(const)
    entries_below_const = round(100.*sum(max_likelihood(:) < const)./136350,2);

% Add a vertical line at 'const'
xline(const, 'r', 'LineWidth', 2, 'Label', sprintf('const = %.2f', const), ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');

% Display the number of entries below 'const'
text(const, max(ylim) * 0.9, sprintf('Major Outlier bin %.2f%%', entries_below_const), ...
    'HorizontalAlignment', 'left', 'Color', 'r', 'FontSize', 10);
end
% Annotate the plot
title('Distribution of Max Participation Probability of all samples');
xlabel('Probability');
ylabel('Probability Frequency');

grid on;
hold off;
    % Define the filename based on k and i
    filename = sprintf([str,'emiter_%d.png'], k);
    
    % Concatenate the save path with the filename
    full_save_path = fullfile(folder_out_path, filename);

    % Save the figure
    saveas(gcf, full_save_path);
    

end
