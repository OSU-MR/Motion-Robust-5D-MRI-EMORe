function [] = fig3d(x, cent_slc, k, i, varargin)

    % Create figure with the given figure number
    figure(varargin{end-1});  % Assuming fignum is always the second to last argument

    % Plot subplots as in the original function
    subplot(3, 4, 1)
    imagesc(extract_slice(abs(x(:,:,:,1,1)), cent_slc(3))); colormap("gray"); axis off;title("Resp 1");
    subplot(3, 4, 2)
    imagesc(extract_slice(abs(x(:,:,:,1,2)), cent_slc(3))); colormap("gray"); axis off;title("Resp 2");
    subplot(3, 4, 3)
    imagesc(extract_slice(abs(x(:,:,:,1,3)), cent_slc(3))); colormap("gray"); axis off;title("Resp 3");
    subplot(3, 4, 4)
    imagesc(extract_slice(abs(x(:,:,:,1,4)), cent_slc(3))); colormap("gray"); axis off;title("Resp 4");
    subplot(3, 4, 5)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,1)), [3 1 2]), cent_slc(2))); colormap("gray"); axis off;
    subplot(3, 4, 6)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,2)), [3 1 2]), cent_slc(2))); colormap("gray"); axis off;
    subplot(3, 4, 7)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,3)), [3 1 2]), cent_slc(2))); colormap("gray"); axis off;
    subplot(3, 4, 8)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,4)), [3 1 2]), cent_slc(2))); colormap("gray"); axis off;
    subplot(3, 4, 9)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,1)), [3 2 1]), cent_slc(1))); colormap("gray"); axis off;
    subplot(3, 4, 10)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,2)), [3 2 1]), cent_slc(1))); colormap("gray"); axis off;
    subplot(3, 4, 11)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,3)), [3 2 1]), cent_slc(1))); colormap("gray"); axis off;
    subplot(3, 4, 12)
    imagesc(extract_slice(permute(abs(x(:,:,:,1,4)), [3 2 1]), cent_slc(1))); colormap("gray"); axis off;

    % Handle optional arguments using varargin
    if nargin > 4 && ~isempty(varargin{1}) && ~isempty(varargin{2})
        nmse = varargin{1};
        ssim = varargin{2};
        sgtitle_str = sprintf("EM Iter %d CS Iter %d NMSE: %.4f SSIM: %.4f", k, i, nmse, ssim);
    elseif nargin > 4 && ~isempty(varargin{1})
        nmse = varargin{1};
        sgtitle_str = sprintf("EM Iter %d CS Iter %d NMSE: %.4f", k, i, nmse);
    elseif nargin > 4 && ~isempty(varargin{2})
        ssim = varargin{2};
        sgtitle_str = sprintf("EM Iter %d CS Iter %d SSIM: %.4f", k, i, ssim);        
    else
        sgtitle_str = sprintf("EM Iter %d CS Iter %d", k, i);
    end

    % Add a supertitle to the figure
    sgtitle(sgtitle_str);

    % Define the filename based on k and i
    filename = sprintf('emiter_%d_cs_iter_%d.png', k, i);

    % Concatenate the save path with the filename (assuming save_path is the last argument)
    full_save_path = fullfile(varargin{end}, filename);

    % Save the figure
    saveas(gcf, full_save_path);

    % Alternatively, you can use exportgraphics for higher quality
    % exportgraphics(gcf, full_save_path, 'Resolution', 300);

end