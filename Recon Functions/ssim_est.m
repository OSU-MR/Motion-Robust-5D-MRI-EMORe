function [avg_ssim,k] = ssim_est(x, x_true)
    % Ensure inputs are the same size
    k = (x_true(:)'*x(:))./norm(x(:)).^2;
    x = k.*x;
    if ~isequal(size(x), size(x_true))
        error('x and x_true must have the same size');
    end
    
    % Initialize variables
    dims = size(x);
    total_ssim = 0;
    count = 0;

    % Compute SSIM slice by slice along the 5th dimension
    for i = 1:dims(5)
        for j = 1:dims(4)
            % Extract slices along 4th and 5th dimensions
            x_slice = squeeze(x(:, :, :, j, i));
            x_true_slice = squeeze(x_true(:, :, :, j, i));
            
            % Compute SSIM for the current slice
            ssim_value = ssim(abs(x_slice./max(abs(x_slice(:)))), abs(x_true_slice./max(abs(x_true_slice(:)))));
            
            % Accumulate SSIM values
            total_ssim = total_ssim + ssim_value;
            count = count + 1;
        end
    end

    % Calculate average SSIM
    avg_ssim = gather(total_ssim / count);
end