function brier_score=wcmp(w_true,w,k,fignum,save_path)
w = squeeze(w);
w_true = squeeze(w_true);
brier_score = mean((w_true - w).^2, 'all');

% Initialize confusion matrix (80 x 80)
confusion_matrix = zeros(size(w_true,2), size(w_true,2));

for n = 1:size(w, 1)
    [~, true_class] = max(w_true(n, :));
    % Add predicted probabilities to the corresponding true class row
    confusion_matrix(true_class, :) = confusion_matrix(true_class, :) + w(n, :);
end

% Normalize by the number of samples in each true class
num_samples_per_class = sum(w_true(:,1:80), 1);
% confusion_matrix = confusion_matrix ./ num_samples_per_class';

% Plot the confusion matrix
figure(fignum);
imagesc(confusion_matrix,[0 max(num_samples_per_class)]);
colorbar;
xlabel('Predicted Class');
ylabel('True Class');
title("EM Iter " + k + " Confusion Matrix"+"  Brier Score: "+ brier_score);

% Define the filename based on k and i
filename = sprintf('weight_emiter_%d.png', k);

% Concatenate the save path with the filename
full_save_path = fullfile(save_path, filename);

% Save the figure
saveas(gcf, full_save_path);
% Brier Score Calculation
end