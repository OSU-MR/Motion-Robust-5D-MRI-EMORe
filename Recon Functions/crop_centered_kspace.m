function kspc_windowed = crop_centered_kspace(kspc_ref, target_shape, tukey_alpha)
    % kspc_ref: Original k-space data
    % target_shape: Cell array specifying desired dimensions, e.g., {96, 96, [], [], []}
    % tukey_alpha: Parameter controlling the tapering of the Tukey window (0 to 1)
    
    % Get the size of the original k-space data
    current_shape = size(kspc_ref);
    num_dims = length(current_shape);
    slices = cell(1, num_dims);
    window = 1; % Initialize the window as 1 (will be updated for each dimension)
    
    for dim = 1:num_dims
        if isempty(target_shape) || length(target_shape) < dim || isempty(target_shape{dim}) || target_shape{dim} >= current_shape(dim)
            % If target size is empty or larger than current size, don't crop this dimension
            slices{dim} = 1:current_shape(dim);
            % Window remains 1 for this dimension
        else
            size_dim = current_shape(dim);
            desired_size = target_shape{dim};
            crop_total = size_dim - desired_size;
            crop_before = floor(crop_total / 2);
            crop_after = crop_total - crop_before;
            start_idx = crop_before + 1;
            end_idx = size_dim - crop_after;
            slices{dim} = start_idx:end_idx;
            
            % Create Tukey window for this dimension
            N = desired_size;
            tukey_win = tukeywin(N, tukey_alpha);
            % Reshape the window to match the dimension
            win_shape = ones(1, num_dims);
            win_shape(dim) = N;
            tukey_win = reshape(tukey_win, win_shape);
            % Multiply the window with the existing window (element-wise)
            window = window .* tukey_win;
        end
    end
    
    % Crop the k-space data
    kspc_cropped = kspc_ref(slices{:});
    
    % Apply the window to the cropped k-space data
    kspc_windowed = kspc_cropped .* window;
end
