function [] = gif3d(x, cent_slc, k, varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gif3d creates a GIF from a 3D (or 5D) dataset 'x' by slicing it and
    % overlaying a super-title containing iteration number, and optionally,
    % PSNR and SSIM values. Additionally, title_str, min_scale, scale_img,
    % and max_scale can be supplied as optional parameters.
    %
    % Usage:
    %   gif3d(x, cent_slc, k, [psnr], [ssim], [title_str], [min_scale],  [max_scale],[scale_img],[sharp_val], save_path, fignum)
    %
    % The last two arguments in varargin must be:
    %   - save_path (a string specifying the folder to save the GIF)
    %   - fignum (a figure number to use)
    %
    % All other optional inputs (if provided) are read in the order given above.
    % If any are omitted or empty, default values are used.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Total number of optional arguments provided
    Nopt_total = numel(varargin);
    if Nopt_total < 2
        error('You must supply at least a save_path and a figure number as the last two varargin elements.');
    end

    % Extract the last two arguments:
    % (Assumed order: second-to-last = save_path, last = fignum)
    save_path = varargin{end-1};
    fignum = varargin{end};

    % Number of optional parameters preceding these two:
    nOpt = Nopt_total - 2;

    % Default values:
    default_title_str = '';
    default_min_scale = min(abs(x(:)));
    default_scale_img = 1;
    default_max_scale = default_scale_img * max(abs(x(:)));

    % Assign optional parameters if provided:
    if nOpt >= 1 && ~isempty(varargin{1})
        psnr_val = varargin{1};
    else
        psnr_val = [];
    end

    if nOpt >= 2 && ~isempty(varargin{2})
        ssim_val = varargin{2};
    else
        ssim_val = [];
    end

    if nOpt >= 3 && ~isempty(varargin{3})
        title_str = varargin{3};
    else
        title_str = default_title_str;
    end

    if nOpt >= 4 && ~isempty(varargin{4})
        min_scale = varargin{4};
    else
        min_scale = default_min_scale;
    end


    if nOpt >= 6 && ~isempty(varargin{6})
        scale_img = varargin{6};
    else
        scale_img = default_scale_img;
    end

    if nOpt >= 5 && ~isempty(varargin{5})
        max_scale = varargin{5};
    else
        max_scale = scale_img * max(abs(x(:)));
    end

    if nOpt >= 7 && ~isempty(varargin{7})
        sharp_val = varargin{7};
    else
        sharp_val = [];
    end
    % Define the filename for the GIF based on k
    filename = title_str+"_iter_"+k+".gif";

    % Build the full save path for the GIF
    full_save_path = fullfile(save_path, filename);

    % Get the size of the time dimension (assumed to be the 4th dimension)
    [~, ~, ~, time_dim, ~] = size(x);

    % Build the super-title base string
    if ~isempty(psnr_val) && ~isempty(ssim_val) && ~isempty(sharp_val)
        sgtitle_base = sprintf("%s Iter %d PSNR: %.4f SSIM: %.4f Sharp: %.4f Frame: ", title_str, k, psnr_val, ssim_val,sharp_val);
    elseif ~isempty(psnr_val)
        sgtitle_base = sprintf("%s Iter %d PSNR: %.4f Frame: ", title_str, k, psnr_val);
    elseif ~isempty(ssim_val)
        sgtitle_base = sprintf("%s Iter %d SSIM: %.4f Frame: ", title_str, k, ssim_val);
    elseif ~isempty(sharp_val)
        sgtitle_base = sprintf("%s Iter %d Sharpness: %.4f Frame: ", title_str, k, sharp_val);
    else
        sgtitle_base = sprintf("%s Iter %d Frame: ", title_str, k);
    end

    % Loop over the time dimension to create the GIF
    for t = 1:time_dim
        % Create figure with the specified figure number
        figure(fignum); 
        clf;

        % Plot subplots for the current time step using the provided scaling
        subplot(3, 4, 1)
        imagesc(extract_slice(abs(x(:,:,:,t,1)), cent_slc(3)), [min_scale max_scale]); 
        colormap("gray"); axis off; title("Resp 1");
        subplot(3, 4, 2)
        imagesc(extract_slice(abs(x(:,:,:,t,2)), cent_slc(3)), [min_scale max_scale]); 
        colormap("gray"); axis off; title("Resp 2");
        subplot(3, 4, 3)
        imagesc(extract_slice(abs(x(:,:,:,t,3)), cent_slc(3)), [min_scale max_scale]); 
        colormap("gray"); axis off; title("Resp 3");
        subplot(3, 4, 4)
        imagesc(extract_slice(abs(x(:,:,:,t,4)), cent_slc(3)), [min_scale max_scale]); 
        colormap("gray"); axis off; title("Resp 4");
        subplot(3, 4, 5)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,1)), [1 3 2]), cent_slc(2)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 6)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,2)), [1 3 2]), cent_slc(2)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 7)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,3)), [1 3 2]), cent_slc(2)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 8)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,4)), [1 3 2]), cent_slc(2)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 9)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,1)), [2 3 1]), cent_slc(1)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 10)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,2)), [2 3 1]), cent_slc(1)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 11)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,3)), [2 3 1]), cent_slc(1)), [min_scale max_scale]); 
        colormap("gray"); axis off;
        subplot(3, 4, 12)
        imagesc(extract_slice(permute(abs(x(:,:,:,t,4)), [2 3 1]), cent_slc(1)), [min_scale max_scale]); 
        colormap("gray"); axis off;

        % Add a super-title to the figure, including the time step (t)
        sgtitle_str = sgtitle_base + string(t);
        sgtitle(sgtitle_str);

        % Capture the frame and write to the GIF file
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if t == 1
            imwrite(imind, cm, full_save_path, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, full_save_path, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
