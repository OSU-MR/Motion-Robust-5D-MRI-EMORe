function [y,k] = psnr_db(ref,x)
    % psnr_db computes the peak signal-to-noise ratio (PSNR) in decibels
    % between a reference signal (or image) 'ref' and a reconstructed signal 'x'
    % after scaling x optimally.
    %
    % Usage:
    %   [y,k] = psnr_db(ref,x)
    %
    % It computes the scaling factor:
    %   k = (ref' * x) / norm(x)^2;
    %
    % and then computes:
    %   y = 20*log10( norm(ref) / norm(ref - k*x) );
    %
    % A higher value of y indicates a better reconstruction.
    
    % Compute the optimal scaling factor
    k = (ref' * x) ./ (norm(x)^2);
    
    % Compute the PSNR in dB
    y = 20 * log10( sqrt(numel(ref)).*max(abs(ref(:))) / norm(ref(:) - k.*x(:)) );
end
