function [w_init, total_unique_pairs, accel] = initWeighting(PEInd, cbin_est, rbin_est, total_samples, opt, mat_size)
%INITWEIGHTING   Initialize EMORe weighting matrices and compute acceleration
%
%   [w_init, total_unique_pairs, accel] = initWeighting(...)
%   Inputs:
%     PEInd         – [N×2] phase/partition indices after SG removal
%     cbin_est      – [N×1] cardiac bin assignments (NaN if rejected)
%     rbin_est      – [N×1] respiratory bin assignments (NaN if rejected)
%     total_samples – scalar, number of remaining readouts
%     opt           – options struct, must contain nCPhases, nRPhases
%     mat_size      – vector of data dims, mat_size(2)=PE1, mat_size(3)=PE2
%
%   Outputs:
%     w_init             – [1×N×1×(nC*nR+1)] initial weight array for EM
%     total_unique_pairs – total count of unique (phase,partition) pairs used
%     accel              – acceleration rate = (PE1*PE2*nC*nR)/total_unique_pairs

% 1) Initialize binary weighting arrays
w_init = zeros(total_samples, opt.nCPhases, opt.nRPhases); % participation weights
wo_init = zeros(total_samples, 1);                         % outlier weights

for i = 1:total_samples
    if isnan(cbin_est(i)) || isnan(rbin_est(i))
        wo_init(i) = 1;  % mark as outlier
    else
        w_init(i, cbin_est(i), rbin_est(i)) = 1;
    end
end
% 2) Compute total sum of unique PE pairs across all bins
total_unique_pairs = 0;
numBins = opt.nRPhases * opt.nCPhases;
for i = 1:numBins
    bin_i = PEInd(w_init(:,i) == 1, :);
    unique_pairs_i = size(unique(bin_i, 'rows'), 1);
    total_unique_pairs = total_unique_pairs + unique_pairs_i;
end
% 3) Display acceleration metric
disp(['Total sum of unique pairs: ', num2str(total_unique_pairs)]);
accel = (mat_size(2)*mat_size(3)*opt.nCPhases*opt.nRPhases) / total_unique_pairs;
disp(['Acceleration rate: ', num2str(accel)]);

% 4) Reshape weights for EM input: [1×N×1×(nBins+1)]
w_init = w_init(:,:);
w_init = [w_init, wo_init];
w_init = reshape(w_init, [1 size(w_init,1) 1 size(w_init,2)]);

disp('Weighting vectors initialized...');
end
