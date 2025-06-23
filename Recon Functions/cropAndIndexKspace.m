function [kspc_ds, PE_lin_ind, y_ind, mat_size] = cropAndIndexKspace(kspc_ds, PEInd, mat_size, opt)
%CROPANDINDEXKSPACE   Crop readouts, combine coils, and compute k-space indices
%
%   [kspc_ds, PE_lin_ind, y_ind, mat_size] = cropAndIndexKspace(
%        kspc_ds, PEInd, mat_size, opt)
%   Inputs:
%     kspc_ds  – 3D k-space data [RO × PE × N] (single-precision array)
%     PEInd    – [N×2] array of phase/partition indices for each readout
%     mat_size – 1×6 vector of data dimensions before cropping
%                 (mat_size(1)=RO, mat_size(2)=PE1, mat_size(3)=PE2, ...)
%     opt      – options struct with fields:
%                   .crop  – [topFrac, bottomFrac] readout cropping fractions
%                   .nc    – number of coils after combination
%
%   Outputs:
%     kspc_ds     – cropped k-space data after coil combining [RO'×PE×N']
%     PE_lin_ind  – GPU array of linear PE indices [N'×1]
%     y_ind       – GPU array of linear k-space subscripts [ (RO'*PE*NC)×1 ]
%     mat_size    – updated dimensions vector after cropping and coil compression

    %% Cropping Readout
    % Determine number of samples to drop from top/bottom of readout
    crop1 = round(opt.crop(1) * mat_size(1));
    crop2 = round(opt.crop(2) * mat_size(1));
    % Find asymetry region size (unused samples)
    asym_size = max(find(kspc_ds(1:end-1,1,1) == 0));
    asym_percent = asym_size / size(kspc_ds,1);
    % Update readout dimension
    mat_size(1) = mat_size(1) - (crop1 + crop2);
    % Perform FFT-based cropping in readout dimension
    kspc_ds = ifftshift(ifft(fftshift(kspc_ds), [], 1));
    kspc_ds = kspc_ds(crop1+1:end-crop2, :, :);
    kspc_ds = fftshift(fft(ifftshift(kspc_ds), [], 1));

    %% Coil combining and indexing
    % Combine coil channels using Walsh/other method
    [kspc_ds, noise_ch] = coilCombine(kspc_ds, opt.nc);
    % Update coil dimension
    mat_size(4) = opt.nc;

    % Convert phase/partition pairs to linear GPU indices
    PE_lin_ind = gpuArray(sub2ind([mat_size(2), mat_size(3)], PEInd(:,1), PEInd(:,2)));

    % Compute linear k-space subscripts for accumulation
    [idy1, idy2, idy3] = ndgrid(1:mat_size(1), PE_lin_ind, 1:mat_size(4));
    y_ind = gpuArray(sub2ind([mat_size(1), mat_size(2)*mat_size(3), mat_size(4)], ...
                   idy1(:), idy2(:), idy3(:)));

    % Clean up temporary variables
    clear idx1 idx2 idx3 idx4 idy1 idy2 idy3 noise_ch;
end
