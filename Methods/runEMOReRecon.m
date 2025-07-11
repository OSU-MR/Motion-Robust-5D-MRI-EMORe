function [x_em, w_em, cost_em, percentage_energy_class] = runEMOReRecon(y, x0, sens, w_init, theta, ...
    y_ind, PE_lin_ind, dims, n_dims, mat_size,cent_slc, ...
    folder_opt_path, folder_gif_path, folder_out_path, timeFull, resp_interp, card_interp, opt)
%RUNEMORERECON   Initialize and execute EMORe 5D MRI reconstruction
%
%   This function performs a 5D free-running cardiac and respiratory motion-resolved
%   MRI reconstruction using the Expectation-Maximization-based Outlier Rejection (EMORe)
%   framework. It takes k-space data and initial parameters to produce a
%   motion-compensated image series.
% --------------------------------------------------------------------------
% It is based on the research article:
% EMORe: Motion‐Robust 5D MRI Reconstruction via
%         Expectation‐Maximization–Guided Binning Correction
%         and Outlier Rejection
%
% TMI Preprint: <arXiv link placeholder>
% Published Abstract: https://doi.org/10.1016/j.jocmr.2024.101509
%
% Author: Syed M. Arshad    | Rizwan Ahmad
% Email: arshad.32@osu.edu  | ahmad.46@osu.edu
%---------------------------------------------------------------------------
%   [x_em, d_em, b_em, w_em, cost_em, cost_em_ax, cost_em_tv, percentage_energy_class] = ...
%       runEMOReRecon(y, x0, sens, w_init, theta, At, a, y_ind, PE_lin_ind,
%                     dims, n_dims, mat_size, cent_slc,
%                     folder_opt_path, folder_gif_path, folder_out_path,
%                     timeFull, resp_interp, card_interp,
%                     fname, opt)
%
%   Inputs:
%     y                – GPU array of k-space data [RO, Readouts, Coils].
%     x0               – Initial image estimate.
%     sens             – Coil sensitivity maps.
%     w_init           – Initial participation weights based on self-gating.
%     theta               – Prior probabilities for bin assignments (theta).
%     y_ind, PE_lin_ind– Linear indexing arrays for k-space sampling.
%     dims, n_dims     – Dimension vectors for Total Variation (TV) operations.
%     mat_size         – Data dimensions [RO, PE1, PE2, Coils, ...].
%     cent_slc         – Center slice index for generating GIF display.
%     folder_opt_path  – Path to save EM optimization logs and parameters.
%     folder_gif_path  – Path to save initialization and final reconstruction GIFs.
%     folder_out_path  – Path to save outlier distribution plots.
%     timeFull         – Full-readout time vector.
%     resp_interp      – Interpolated respiratory trace.
%     card_interp      – Interpolated cardiac trace.
%     opt              – Struct with optimization options:
%                        .scale:     Image scaling factor.
%                        .sig_scl:   Noise standard deviation (sigma).
%                        .tau:       Outlier rejection threshold.
%                        .gStp:      Gradient descent step size for M-step.
%                        .mu:        ADMM penalty parameter.
%                        .lam:       TV regularization hyperparameters.
%                        .iIter1:    Number of ADMM iterations for initialization.
%                        .iIter2:    Frequency of E-step updates (in M-step iterations).
%                        .vrb:       Verbosity, iteration frequency to display progress.
%                        .thresh:    Convergence threshold for image difference.
%                        .emIter:     Maximum number of EMORe iterations.
%                        .scale_img: Image scaling for GIF visualization.
%
%   Outputs:
%     x_em, d_em, b_em, w_em – Final EMORe variables (image, TV aux vars, weights).
%     cost_em, cost_em_ax, cost_em_tv – Cost history vectors (total, data fidelity, TV).
%     percentage_energy_class       – Percentage of readouts assigned to each respiratory state and the outlier bin.

% --- INITIALIZATION OF COST VECTORS ---
% Pre-allocate arrays to store cost function values at each iteration for monitoring convergence.
cost_em = zeros(opt.emIter,1);       % Total cost (Data Fidelity + TV)
cost_em_ax = zeros(opt.emIter,1);    % Data fidelity component of the cost
cost_em_tv = zeros(opt.emIter,1);    % Total variation (TV) component of the cost

% --- INITIALIZATION OF RECONSTRUCTION VARIABLES ---
% Initialize variables for the reconstruction, moving them to the GPU for acceleration.
w_sqrt = gpuArray(single(sqrt(w_init))); % Initial weights (sqrt is used as forward operator 'a' applies it)
x = gpuArray(repmat(x0,[1 1 1 20 4])); % Replicate initial image for all cardiac/respiratory bins
gradA = zeros(size(x),'like',y);      % Gradient of the data fidelity term
x = x./opt.scale;                     % Scale initial image estimate
L = size(y,1)*size(y,3);                         % L: Collective length of a readout from all coils
sig_sq = opt.sig_scl.^2;               % σ^2 Noise variance 
% Define forward (A) and adjoint (At) operators as function handles for convenience.
At = @(y,w_sqrt) At_operator(y,mat_size(1:4),sens,y_ind,w_sqrt);
A  = @(x,w_sqrt) A_operator(x,[mat_size(1:4),1,1],sens,PE_lin_ind,w_sqrt);

% Initialize auxiliary variables 'd' and 'b' for the ADMM solver for TV minimization.
d = zeros([dims,n_dims],'like',y); % Split-Bregman variable for TV
b = zeros([dims,n_dims],'like',y); % Bregman variable for TV

% --- IMAGE INITIALIZATION (PARTIAL M-STEP) ---
% This section corresponds to the initialization described in Section II-A-3 of the paper.
% It partially solves the M-step optimization problem (Eq. 2) using the initial
% SG-based weights (Eq. 3) to get a reasonable starting image x^(0).
disp("EMORe Image Initialization Started....");
stop = 0; % Loop control flag
i = 0;    % Iteration counter
while (~stop)
        i = i+1;
        % (ADMM): performs gradient descent on the image 'x'.
        % The choice of 4 iterations is a practical one for the ADMM subproblem.
        for c = 1:4 
            % Calculate data fidelity gradient for all motion bins.
            for j=1:size(x(:,:,:,:),4)
                gradA(:,:,:,j) =  (1./(L.*sig_sq)).*At(A(x(:,:,:,j),w_sqrt(:,:,:,j))- (w_sqrt(:,:,:,j).*y),w_sqrt(:,:,:,j));
            end
            % Calculate the gradient of the TV regularization term.
            gradT = opt.mu.* adjTV(TV(x,dims)-d+b,dims);
            % Update the image 'x' by taking a step in the negative gradient direction.
            x = x - opt.gStp.*(gradA+gradT);
        end
        
        % ADMM updates for auxiliary variables 'd' and 'b' (proximal operator for TV).
        tvx = TV(x,dims); % TV operator
        d = sth(tvx + b, opt.lam./opt.mu); % Soft-thresholding step
        b = b + (tvx - d); % Bregman variable update
        
        % Check for stopping criterion.
        if(i==opt.iIter1)
            stop = 1;
        end
end
% Save a GIF of the initialized image for visual inspection.
gif3d(x,cent_slc,0,[],[],'EMORe Image Initialization',[],[],opt.scale_img,[],folder_gif_path,14);
disp("EMORe Image Initialization Completed...")

%% --- EMORe RECONSTRUCTION ---
% This is the main loop of the EMORe algorithm, which iteratively performs
% the E-step and M-step until convergence as outlined in Algorithm 1.
disp("-----------------------------------------------------------")
disp("EMORe Reconstruction Started...")
tau_sq = (opt.sig_scl*opt.tau_factor).^2; % Pre-calculate outlier threshold term from Eq. (1c)
percentage_energy_class = zeros(5,1); % Array to store readout energy per bin
w_sqrt = gpuArray(single(sqrt(w_init)));   % Reset weights to initial SG-based assignment
r = zeros(1,size(y,2),1,mat_size(5)*mat_size(6),'like',y); % Residual error array
likelihood = zeros([size(w_sqrt,1),size(w_sqrt,2),size(w_sqrt,3), size(w_sqrt, 4)],'like',w_sqrt);

% Set up a log file to save progress.
logFile = fullfile(folder_opt_path, [opt.fname,'_EM.txt']);
fid = fopen(logFile, 'a');  % Open for appending

% Initialize loop control variables.
stop = 0; % Loop control flag
k = 0;    % M-step iteration counter
e = 0;    % E-step iteration counter
xi = zeros(size(x),'like',x); % Stores image from previous E-step for convergence check
t1 = tic; % Start timer for the entire reconstruction

while (~stop)
    t2 = tic; % Start timer for the current iteration
    k = k+1; % Increment M-step counter

    % --------------------- E-STEP ------------------------------------------
    % This block corresponds to Section II-A-1 and Eq. (1). It updates the
    % participation weights 'w' based on the current image estimate 'x'.
    % The E-step is performed every 'opt.iIter2' M-step iterations.
    if(mod(k,opt.iIter2)==1)
        e = e+1; % Increment E-step counter
        
        % Calculate the squared residual norm ||A(n,k)x - y_n||^2 for each readout.
        for j=1:size(x(:,:,:,:),4)
            r(:,:,:,j) = sum(sum(abs(A(x(:,:,:,j),[])-y).^2,1),3);
        end
        
        % Calculate the likelihood p(y_n|x, r_n=k) for valid bins (Eq. 1b)
        likelihood(:,:,:,1:80) = (exp(-r./(L.*sig_sq)).*theta(:,:,:,1:80))./(pi.*sig_sq);
        % Calculate the likelihood for the outlier bin (Eq. 1c)
        likelihood(:,:,:,81)  =  (exp(-tau_sq/(sig_sq)).*theta(:,:,:,81))./(pi.*sig_sq);
        
        % Update weights for the M-step. The forward model uses sqrt(w).
        % Normalize to get posterior probabilities (the new weights), as in Eq. (1a).
        w_sqrt = sqrt(likelihood./sum(likelihood,4));
                
        % Calculate the normalized difference between images from consecutive E-steps.
        norm_diff_xem = norm(x(:)-xi(:))./norm(xi(:));
        xi = x; % Store current image for the next E-step's convergence check
    end

    % --------------------- M-STEP ------------------------------------------
    % This block corresponds to Section II-A-2 and Eq. (2). It updates the
    % image estimate 'x' using the new weights 'w' from the E-step.
    % This is a partial solve (GEM approach) with 4 inner iterations.
    for c = 1:4 
        % Calculate data fidelity gradient for all motion bins using updated weights.
        for j = 1:size(x(:,:,:,:),4)
            gradA(:,:,:,j) =  (1./(L.*sig_sq)).*At(A(x(:,:,:,j),w_sqrt(:,:,:,j))- (w_sqrt(:,:,:,j).*y),w_sqrt(:,:,:,j));
        end
        % Calculate the TV regularization gradient.
        gradT = opt.mu.* adjTV(TV(x,dims)-d+b,dims);
        % Update the image 'x' via gradient descent.
        x = x - opt.gStp.*(gradA+gradT);
    end
    
    % ADMM updates for auxiliary variables 'd' and 'b'.
    tvx = TV(x,dims);
    d = sth(tvx + b, opt.lam./opt.mu);
    b = b + (tvx - d); 
    
    % --- LOGGING AND CONVERGENCE CHECK ---
    % Periodically display and log progress.
    if(mod(e,opt.vrb)==0 || e==1)
        cost_em_ax(k) = round((1./(L.*sig_sq)).*norm((L.*sig_sq).*gradA(:)).^2);
        cost_em_tv(k) = round(sum(opt.lam.*sum(abs(reshape(tvx,[],n_dims)),1)));
        cost_em(k) = cost_em_ax(k)+cost_em_tv(k);
        % Form and display the log message.
        logMessage = sprintf('EMORe Iter %d Time/iter(s) = %.2f Cost: %d Cost ax: %d Cost tv: %d Diff: %.4f', ...
            e, toc(t2), cost_em(k),cost_em_ax(k),cost_em_tv(k), norm_diff_xem);
        disp(logMessage);
        fprintf(fid, '%s\n', logMessage); % Write to log file
    end
    
    % Check for algorithm convergence.
    if(norm_diff_xem < opt.thresh || e==opt.emIter)
            stop = 1;
    end
end
% --- FINALIZATION ---
% Stop timer and log final results.
t_em = toc(t1)./60; % Total reconstruction time in minutes
logMessage = sprintf('EMORe is completed in %.2f mins %s', t_em,datestr8601([],'*ymdHMS'));
disp(logMessage);
fprintf(fid, '%s\n', logMessage);
fclose(fid);  % Close the log file

% --- SAVE RESULTS ---
% Generate GIF of final reconstructed image series.
gif3d(x,cent_slc,k,[],[],'EMORe',[],[],opt.scale_img,[],folder_gif_path,14);
% Gather data from GPU to CPU.
x_em = gather(x);
w_em = gather(w_sqrt.^2);
% Save final variables and parameters to a .mat file.
save(fullfile(folder_opt_path,[opt.fname,'_iter',num2str(k),'_em.mat']), 'x_em','opt', ...
    'norm_diff_xem','percentage_energy_class','cost_em','t_em','k','-v7.3');

%% --- PLOTTING ---
% Plot physiological traces overlaid with the outlier bin assignments to visualize
% which parts of the cardiac/respiratory cycles were rejected.
figure;
% Plot Respiratory Signal vs. Outlier Weights
subplot(2,1,1);
yyaxis left;
plot(timeFull, resp_interp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Resp Interp');
ylabel('Respiration Signal');
yyaxis right;
bar(timeFull, w_sqrt(1,:,1,81), 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 1);
ylim([0 1]);
ylabel('Outlier Bin Assignment % (w_out)');
title(['Respiration Signal with ' num2str(round(percentage_energy_class(5),2)) '% Rejected Outliers']);
xlabel('Time (s)');
legend("Resp","Outliers");
grid on;

% Plot Cardiac Signal vs. Outlier Weights
subplot(2,1,2);
yyaxis left;
plot(timeFull, card_interp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Card Interp');
ylabel('Cardiac Signal');
yyaxis right;
bar(timeFull, w_sqrt(1,:,1,81), 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 1);
ylim([0 1]);
ylabel('Outlier Bin Assignment % (w_out)');
title(['Cardiac Signal with ' num2str(round(percentage_energy_class(5),2)) '% Rejected Outliers']);
xlabel('Time (s)');
legend("Card","Outliers");
grid on;

% Save the generated figure.
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(folder_out_path, 'outliers.png'));

end