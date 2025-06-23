% Hermitian 3D Operator
function y = multTr(obj,x)

    ndims_x = numel(obj.frame_size);
    y = zeros([obj.frame_size(1:3),size(obj.maps,4),obj.frame_size(4:ndims_x)],'single','gpuArray');

    

    
    % Weigthed Least Squares
    if ~isempty(obj.mask_weights)
        x = obj.mask_weights.*x;
    end
    
    
    
    y(obj.mask_patterns) = x;
    % iFFT
    y = sqrt(size(y,1)*size(y,2)*size(y,3))*ifft(ifft(ifft(y,[],1),[],2),[],3);
    % Coil Multiply
    
    % Check if multiple sensitivity maps are used for espirit
    if(size(obj.maps,5) == 1)
        y = y.*conj(obj.maps);
    else
        y = cat(ndims_x+2,y.*conj(obj.maps(:,:,:,:,1)),y.*conj(obj.maps(:,:,:,:,2)));
    end
      
    y = squeeze(sum(y,4));

end