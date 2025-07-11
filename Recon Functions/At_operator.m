% Hermitian 3D Operator
function x = At_operator(y,mat_size,sen,lin_ind,w)

    if nargin < 5 || isempty(w)
        % 'w' is not given, initialize to ones
        performWStep = false; % Flag to indicate whether to perform the step with 'w'
    else
        performWStep = true; % Flag to indicate that 'w' is given
    end
    if(performWStep)
        y = w.*y;
    end
    x = accumarray(lin_ind,y(:),[prod(mat_size) 1]);
    % if(avg)
    % w_norm = accumarray(PE_lin_ind,w(:).^2,[prod(mat_size(2:3)) 1]);
    % w_norm(w_norm == 0) = inf;
    % x = reshape(x,mat_size)./reshape(w_norm,[1,mat_size(2),mat_size(3)]);
    % else
    x = reshape(x,mat_size);   
    % end
    % iFFT
    x = ifft3_shift(x);  
    % x = ifft3(x);  
    x = x.*conj(sen);
    x = squeeze(sum(x,4));

end