% Forward 3D Operator
function y = A(x,mat_size,maps,PE_lin_ind,w)
    if nargin < 5 || isempty(w)
        % 'w' is not given, initialize to ones
        performWStep = false; % Flag to indicate whether to perform the step with 'w'
    else
        performWStep = true; % Flag to indicate that 'w' is given
    end
    % Image domain multiply and FFT
    ndims_x = ndims(x);
    x = permute(x, [1:3,ndims_x+1,4:ndims_x]);    
    x = x.*maps;
    x = fft3_shift(x);
    % x = fft3(x);


    % Downsample
    x = reshape(x,mat_size(1),mat_size(2)*mat_size(3),mat_size(4),[]);
    y = x(:,PE_lin_ind,:,:,:);
    if(performWStep)
        y = w.*y;
    end
    

end