function [ kdata_coil, noise_ch] = coilCombine( kdata_coil, n_channels)
    % Channel compression =================================================
    n_coils = size(kdata_coil,2);
    if(n_coils == n_channels)
        noise_ch = kdata_coil(:,end,:);
        disp("Coil compression not needed. Using last coil for noise")
    elseif  n_coils > n_channels
        fprintf('Compressing %d channels to %d\n',size(kdata_coil,2),n_channels);
        dims = ndims(kdata_coil);  % Get the number of dimensions in the array
        order = [1,3,2];  % Move coil (4th) dimension to the end
        kdata_coil = permute(kdata_coil, order);  % Rearrange the dimensions
        N = size(kdata_coil);
        % Performing SVD for Channel Compression (optional) 
        kdata_coil = reshape(kdata_coil,[],N(dims));
        [~,S,V] = svd(kdata_coil'*kdata_coil,0);
        % Get the noise channel
        noise_ch = kdata_coil*V(:,:);
        noise_ch = noise_ch(:,end-3:end);
        noise_ch = reshape(noise_ch,N(1),N(2),[]);
        noise_ch = squeeze(noise_ch);

        kdata_coil = kdata_coil*V(:,1:n_channels);
        kdata_coil = reshape(kdata_coil,[N(1:dims-1),n_channels]);
        kdata_coil = squeeze(kdata_coil);

        % Compute the inverse of the permutation
        inverse_order(order) = 1:length(order);  % Inverse of the permutation
        % Rearrange to the original order
        kdata_coil = permute(kdata_coil, inverse_order);
    end


end