function [x_cs, cost_cs] = runCSRecon(y, x0, sens, w_init, y_ind, ...
    PE_lin_ind, dims, n_dims, mat_size,cent_slc, folder_opt_path, folder_gif_path, opt)
%RUNCSRRECON   Execute a standard Compressed Sensing (CS) 5D MRI reconstruction.
%
%   This function performs a 5D free-running cardiac and respiratory motion-resolved
%   MRI reconstruction using a conventional compressed sensing approach with
%   Total Variation (TV) regularization. The binning is fixed based on the
%   initial self-gating signal and is not updated.
%-------------------------------------------------------------------------------
% Author: Syed M. Arshad
% Email: arshad.32@osu.edu
%-------------------------------------------------------------------------------
%   [x_cs, cost_cs] = runCSRecon(y, x0, sens, w_init, y_ind, PE_lin_ind, ...
%                                dims, n_dims, mat_size, cent_slc, ...
%                                folder_opt_path, folder_gif_path, opt)
%
%   Inputs:
%     y                – GPU array of k-space data [RO, Readouts, Coils, Bins].
%     x0               – Initial image estimate.
%     sens             – Coil sensitivity maps.
%     w_init           – Initial, fixed participation weights from self-gating.
%     y_ind, PE_lin_ind– Linear indexing arrays for k-space sampling.
%     dims, n_dims     – Dimension vectors for Total Variation (TV) operations.
%     mat_size         – Data dimensions [RO, PE1, PE2, Coils, ...].
%     cent_slc         – Center slice index for generating GIF display.
%     folder_opt_path  – Path to save optimization logs and parameters.
%     folder_gif_path  – Path to save reconstruction GIF.
%     opt              – Struct with optimization options:
%                        .scale:     Image scaling factor.
%                        .gStp:      Gradient descent step size.
%                        .mu:        ADMM penalty parameter.
%                        .lam:       TV regularization hyperparameters.
%                        .vrb:       Verbosity, iteration frequency to display progress.
%                        .thresh:    Convergence threshold for image difference.
%                        .csIter:     Maximum number of iterations.
%                        .scale_img: Image scaling for GIF visualization.
%
%   Outputs:
%     x_cs             – Final reconstructed image series.
%     cost_cs          – Cost history vector.

    % --- INITIALIZATION OF COST VECTORS ---
    % Pre-allocate arrays to store cost function values for monitoring convergence.
    cost_cs = zeros(opt.csIter,1);       % Total cost (Data Fidelity + TV)
    cost_cs_ax = zeros(opt.csIter,1);    % Data fidelity component of the cost
    cost_cs_tv = zeros(opt.csIter,1);    % Total variation (TV) component of the cost

    % --- INITIALIZATION OF RECONSTRUCTION VARIABLES ---
    % Initialize variables and move them to the GPU for accelerated computation.
    w = gpuArray(single(sqrt(w_init))); % Use fixed initial weights (sqrt is for the forward model)
    x = gpuArray(repmat(x0,[1 1 1 20 4])); % Replicate initial image for all cardiac/respiratory bins
    gradA = zeros(size(x),'like',y);      % Pre-allocate gradient of the data fidelity term
    x = x./opt.scale;                     % Scale initial image estimate

    % Define forward (A) and adjoint (At) operators as function handles.
    at = @(y,w) At(y,mat_size(1:4),sens,y_ind,w);
    a  = @(x,w) A(x,[mat_size(1:4),1,1],sens,PE_lin_ind,w);
    
    % Initialize auxiliary variables 'd' and 'b' for the ADMM solver.
    d = zeros([dims,n_dims],'like',y); % Split-Bregman variable for TV
    b = zeros([dims,n_dims],'like',y); % Bregman variable for TV
    
    disp("CS Algorithm Started....");
    
    % Generate a GIF of the time-averaged initial image.
    gif3d(x(:,:,:,1,:),cent_slc,0,[],[],'Time-averaged',[],[],opt.scale_img,[],folder_gif_path,14);
    
    % --- SETUP LOGGING ---
    % Create and open a log file to record the progress of the reconstruction.
    logFile = fullfile(folder_opt_path, [opt.fname,'_CS.txt']);
    fid = fopen(logFile, 'a');  % open for appending
    
    % --- MAIN RECONSTRUCTION LOOP (ADMM) ---
    % This loop iteratively solves the TV-regularized least-squares problem.
    stop = 0; % Loop control flag
    i = 0;    % Iteration counter
    t1 = tic; % Start timer for the entire reconstruction
    
    while (~stop)
            i = i+1; % Increment iteration counter
            t2 = tic; % Start timer for the current iteration
            xi = x;   % Store current image estimate to check for convergence

            % Inner loop performs gradient descent on the image 'x' (x-minimization step of ADMM).
            % This is a partial solve with a fixed number of sub-iterations (4).
            for c = 1:4 
                % Calculate the gradient of the data fidelity term for all motion bins.
                for j=1:size(x(:,:,:,:),4)
                    gradA(:,:,:,j) =  at(a(x(:,:,:,j),w(:,:,:,j))- (w(:,:,:,j).*y),w(:,:,:,j));
                end
                % Calculate the gradient of the TV regularization term.
                gradT = opt.mu.* adjTV(TV(x,dims)-d+b,dims);
                % Update the image 'x' by taking a step in the negative gradient direction.
                x = x - opt.gStp.*(gradA+gradT);
            end

            % ADMM updates for auxiliary variables 'd' and 'b' (proximal operator for TV).
            tvx = TV(x,dims); % Apply TV operator
            d = sth(tvx + b, opt.lam./opt.mu); % Soft-thresholding step
            b = b + (tvx - d); % Bregman variable update
            
            % Calculate the normalized difference between consecutive image estimates.
            norm_diff_xcs = norm(x(:)-xi(:))./norm(xi(:));
           
            % --- LOGGING AND CONVERGENCE CHECK ---
            % Periodically display and log progress.
            if(mod(i,opt.vrb)==0 || i==1)
                cost_cs_ax(i) = round(0.5.*norm(gradA(:)).^2);
                cost_cs_tv(i) = round(sum(opt.lam.*sum(abs(reshape(tvx,[],n_dims)),1)));
                cost_cs(i) = cost_cs_ax(i)+cost_cs_tv(i);
    
                % Form and display the log message.
                logMessage = sprintf('CS Iter %d Time/iter(s) = %.2f Cost: %d Cost ax: %d Cost tv: %d Diff: %.4f', ...
                    i, toc(t2), cost_cs(i),cost_cs_ax(i),cost_cs_tv(i), norm_diff_xcs);
                disp(logMessage);
                fprintf(fid, '%s\n', logMessage); % Write to log file
            end

            % Check for algorithm convergence based on image difference or max iterations.
            if(norm_diff_xcs < opt.thresh || i == opt.csIter)
                stop = 1;
            end
    end
    % --- FINALIZATION AND SAVING ---
    % Stop timer and log the final results.
    t_cs = toc(t1)./60; % Total reconstruction time in minutes
    
    % Recalculate cost for the final iteration.
    cost_cs_ax(i) = round(0.5.*norm(gradA(:)).^2);
    cost_cs_tv(i) = round(sum(opt.lam.*sum(abs(reshape(tvx,[],n_dims)),1)));
    cost_cs(i) = cost_cs_ax(i)+cost_cs_tv(i);
    
    % Form and log the final message.
    logMessage = sprintf('CS Iter %d Time/iter(s) = %.2f Cost: %d Cost ax: %d Cost tv: %d Diff: %.4f', ...
        i, toc(t2), cost_cs(i),cost_cs_ax(i),cost_cs_tv(i), norm_diff_xcs);
    disp(logMessage);
    fprintf(fid, '%s\n', logMessage);
    logMessage = sprintf('CS is completed in %.2f mins %s', t_cs,datestr8601([],'*ymdHMS'));
    disp(logMessage);
    fprintf(fid, '%s\n', logMessage);    
    fclose(fid);  % Close the log file
    
    % Generate a GIF of the final reconstructed image series.
    gif3d(x,cent_slc,i,[],[],'CS',[],[],opt.scale_img,[],folder_gif_path,14);
    
    % Gather data from GPU to CPU for saving.
    x_cs = gather(x);
    
    % Save final variables and parameters to a .mat file.
    save(fullfile(folder_opt_path, [opt.fname,'_iter',num2str(i),'_cs.mat']), 'x_cs','opt','norm_diff_xcs','t_cs','i','-v7.3');
end