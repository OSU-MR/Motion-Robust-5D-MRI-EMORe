% Hermitian 3D Operator
function x = At2(y, mat_size, sen, ind_loc, w)

    if nargin < 5 || isempty(w)
        % 'w' is not given, initialize to ones
        performWStep = false; % Flag to indicate whether to perform the step with 'w'
    else
        performWStep = true; % Flag to indicate that 'w' is given
    end
    
    if performWStep
        y = w .* y;
    end
    
    y = permute(y, [1, 3, 4, 2]);
    x = zeros(mat_size(1), mat_size(4), mat_size(5) * mat_size(6), mat_size(2) * mat_size(3));
    
    % Use accumarray to replace the inefficient for-loop
    % Flatten y to a 2D matrix for easier indexing with accumarray
    y_flat = reshape(y, [], size(ind_loc, 1));
    
    % Create linear indices for accumarray
    linear_indices = repmat(ind_loc(:)', size(y_flat, 1), 1);
    
    % Accumulate values from y_flat into x_flat
    x_flat = accumarray(linear_indices(:), y_flat(:), [size(y_flat,1),mat_size(2) * mat_size(3)], @sum);
    
    % Reshape x_flat back to the original dimensions
    x = reshape(x_flat, [mat_size(1), mat_size(4), mat_size(5) * mat_size(6), mat_size(2) * mat_size(3)]);
    
    % Normalize by the counter (using accumarray to count occurrences)
    counter = accumarray(ind_loc, 1, [mat_size(2) * mat_size(3), 1]);
    counter(counter == 0) = inf;
    x = x ./ permute(counter, [2 3 4 1]);
    
    % Reshape and permute back to the original form
    x = permute(reshape(x, [mat_size(1), mat_size(4), mat_size(5), mat_size(6), mat_size(2), mat_size(3)]), [1, 5, 6, 2, 3, 4]);

    % iFFT
    x = ifft3_shift(x);  
    x = x .* conj(sen);
    x = squeeze(sum(x, 4));

end
