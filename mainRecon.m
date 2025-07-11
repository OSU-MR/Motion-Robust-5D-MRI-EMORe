%==========================================================================
% EMORe 5D MRI Reconstruction Main Script
%
% EMORe: Motion‐Robust 5D MRI Reconstruction via
%         Expectation‐Maximization–Guided Binning Correction
%         and Outlier Rejection
%
% TMI Preprint: <arXiv link placeholder>
% Published Abstract: https://doi.org/10.1016/j.jocmr.2024.101509
%
% Authors: Syed M. Arshad, Lee C. Potter, Chong Chen,
%          Preethi C. Subramanian, Yingmin Liu,
%          Christopher Crabtree, Xuan Lei, Juliet Varghese,
%          Matthew S. Tong, Rizwan Ahmad
% Email: arshad.32@osu.edu
%
% Brief:
%   1) Prompts user to select a raw .mat dataset from ./Datasets/
%   2) Loads self‐gated 5D k‐space and metadata
%   3) Sets reconstruction parameters for CS and EMORe
%   4) Hands off to the EMORe pipeline (in Recon Functions/)
%==========================================================================

%% --- Setup and Paths ----------------------------------------------------
addpath(genpath('./Recon Functions/'));   % include all helper functions
addpath(genpath('./Methods/'));   % include all helper functions

close all;
clear;

%% --- Import the Data ---------------------------------------------------
% Prompt user to choose raw data .mat files in ./Datasets/
[fileName, filePath] = uigetfile({'*.mat'}, ...
    'Select raw data file', './Datasets');
if isequal(fileName,0)
    error('No file selected. Exiting.');
end
fullFileName = fullfile(filePath, fileName);

disp("----------------------------------------------------------------");
disp("  EMORe Reconstruction");
disp("  Loading raw data: " + fileName);
disp("----------------------------------------------------------------");
load(fullFileName);   % must load variables: Lin, Par, sqzSize, kspc_ds, samp, param
disp("  Raw data import successful.");
disp("----------------------------------------------------------------");

%% --- Reconstruction Parameters -----------------------------------------
disp("  Setting reconstruction parameters (see manuscript for details)...");

% Comparison flag
opt.cs         = 1;                   % 1 = run CS recon; 0 = skip CS

% ADMM tuning
opt.lam        = 2e-2 .* [1 1 1 5 3];  % [λ_sx, λ_sy, λ_sz, λ_c, λ_r]
opt.gStp       = 1e-3;                % gradient step size (γ)
opt.mu         = 20;                % Lagrange multiplier (μ)

% EMORe priors
opt.out_factor = 0.05;                % α_o (outlier prior)
opt.sg_factor  = 0.85;                % α_g (self‐gating prior)

% Iteration counts
opt.emIter      = 60;                 % Max EMORe iterations (J)
opt.iIter1     = 8;                   % inner ADMM iterations I₁
opt.iIter2     = 4;                   % inner ADMM iterations I₂

opt.csIter      = 240;                 % Max CS iterations (J)


% Outlier threshold & noise level
opt.tau_factor        = 3;                   % τ-factor (τ/σ) (outlier threshold)
opt.sig               = 1e-5;                % σ (noise std; data pre‐scaled)

% Convergence
opt.thresh     = 1e-4;                % η (convergence criterion)

% Visualization of results
opt.scale_img  = 0.7;                 % clipping factor for output images/GIFs
opt.vrb = 4;                          % iterations between progress (cost) updates

% Binning
opt.nRPhases   = 4;                   % respiratory bins
opt.nCPhases   = 20;                  % cardiac bins

% Coil compression & filtering
opt.nc         = 8;                   % # coils post‐compression
opt.rBPF       = [0.02, 0.5];         % resp. bandpass [Hz]
opt.cBPF       = [0.5,  3.5];         % card. bandpass [Hz]

% Readout cropping (% from top, bottom) e.g. [0,0]=no crop; [0.25,0.25]=50% total
opt.crop       = [0.40, 0.35];        % used for Volunteer#1

disp("  Parameters configured:");
disp(opt);
disp("----------------------------------------------------------------");


%% --- Setup output directories for reconstruction -----------------------
% tag            – optional suffix to distinguish runs (e.g. '_test1')
% fname          – base name of the dataset (without “.mat”)
% fname_datetime – unique folder name combining fname + timestamp

tag = '';  
opt.fname = fileName(1:end-4);                               % strip “.mat”
timestamp = datestr8601([], '*ymd', 'HM');
datetimeDir    = fullfile('recons', opt.fname, [opt.fname '_' timestamp tag]);
fname_datetime = [opt.fname '_' timestamp tag];

% subfolders for gifs, outlier plots, and parameter logs
folder_gif  = fullfile(datetimeDir, 'gif');
folder_out  = fullfile(datetimeDir, 'outliers');
folder_opt  = fullfile(datetimeDir, 'params');

% create the directory tree in one go
mkdir(datetimeDir);
mkdir(folder_gif);
mkdir(folder_out);
mkdir(folder_opt);

% store full paths for later use
folder_parent_path = fullfile(pwd, 'recons', opt.fname);
folder_study_path  = datetimeDir;
folder_gif_path    = folder_gif;
folder_out_path    = folder_out;
folder_opt_path    = folder_opt;

disp("New folders for reconstructions created:");
disp("  " + folder_study_path);
disp("  - GIFs:     " + folder_gif_path);
disp("  - Outliers: " + folder_out_path);
disp("  - Params:   " + folder_opt_path);
disp("----------------------------------------------------------------");

%% GPU Selection

% Script to select a GPU interactively, displaying information similar to nvidia-smi.
% Only prompts for selection if multiple GPUs are available.
[selectedGPUIndex, gpuInfo] = selectGPUInteractive();
% NOTE: This reconstruction code is optimized for GPU execution.
%       If you must run on CPU only, make appropriate changes accordingly


%% --- Acquisition Indexing & Self-Gating Initialization -------------
% Builds k-space readout indices, identifies the self-gating channel,
% crops pre-cycling data, and computes data dimensions.

% 1) Phase/partition pairs for each readout
PEInd_SG(:,1) = Lin; % Phase encoding indices
PEInd_SG(:,2) = Par; % Partition encoding indices

% 2) Identify the dominant SG pair and sampling interval
% this function automatically detects SG signal in raw data and its frequency
% Note: if your data does not contain SG or you have pre-binned data, bypass these steps
[opt.SG_pair, opt.sgI, opt.SGfirstOccurrence, opt.SGlastOccurrence] = ...
    estimateSelfGating(PEInd_SG);

% 3) Cropping data to have complete SG cycles from first to last SG sample
PEInd_SG = PEInd_SG(opt.SGfirstOccurrence+1:opt.SGlastOccurrence,:); 
kspc_ds = kspc_ds(:,:,opt.SGfirstOccurrence+1:opt.SGlastOccurrence);
total_samples = size(PEInd_SG,1); % Total readout samples
SG_locs_ind = opt.sgI:opt.sgI:total_samples; % Location indices for SG Readouts
total_SG_samples = numel(SG_locs_ind); % Total SG samples 

% mat_size holds the dimensions of our 5D data volume:
%   dim 1 = readout samples
%   dim 2 = phase‐encode steps
%   dim 3 = partition‐encode steps
%   dim 4 = coil channels
%   dim 5 = cardiac bins
%   dim 6 = respiratory bins

% 4) Store 5D volume dimensions
mat_size = [ ...
    size(kspc_ds,1), ...               % readout length
    sqzSize(3), sqzSize(4), ...        % phase & partition
    sqzSize(2), ...                    % coils
    opt.nCPhases, opt.nRPhases         % cardiac & respiratory bins
];

disp("  Data dimensions: [" + sprintf('%d×%d×%d×%d×%d', ...
    mat_size(1), mat_size(2), mat_size(3), mat_size(4), mat_size(5)) + "]");

disp("Acquisition and SG setup complete:");
disp("  Total readouts:       " + total_samples);
disp("  SG interval (sgI):    " + opt.sgI);
disp("  SG samples:           " + total_SG_samples);
disp("  Data mat_size:        [" + sprintf('%d×%d×%d×%d×%d×%d', mat_size) + "]");
disp("----------------------------------------------------------------");

%% Pre-processing--SG Signal Extraction and Initial SG-based binning

% 1) Extract self-gating
SG = fftshift(ifft(ifftshift(kspc_ds(:,:,SG_locs_ind), 1), [], 1), 1);
opt.TR = param.TR./1e6; % Converting the units of Repetition Time (TR) to seconds
opt.TR_SG = opt.TR.*opt.sgI;  % Calculating Repetition time for SG signal
opt.FS_SG = 1/opt.TR_SG; % SG signal frequency
timeFull = (0:1:total_samples-1)*(opt.TR); % Time instances for readout acqusition starting at 0
SG_timeFull = timeFull(SG_locs_ind); % Time instances for SG acqusition starting at 0


% 2) Build Casorati and filter into resp/card signals
[resp, card] = filterSGSignals(SG, opt, folder_opt_path);

% Plot Card and Resp Signals
figure;
subplot(211); 
% Plot `resp_true` with a solid blue line and `resp` with a dashed red line
plot(SG_timeFull, card, 'r-', 'LineWidth', 1); 
% Add labels and legend for better clarity
xlabel('Time');
ylabel('Amplitude');
legend('SG Extracted', 'Location', 'best');
title('Cardiac Signal');
hold off;


subplot(212); 
% Plot `resp_true` with a solid blue line and `resp` with a dashed red line
plot(SG_timeFull, resp, 'r-', 'LineWidth', 1); 
% Adjust the horizontal lines in `bin_lines` with a dotted green line

% Add labels and legend for better clarity
xlabel('Time');
ylabel('Amplitude');
legend( 'SG Extracted', 'Location', 'best');
title('Resp Respiratory Signal');
hold off;


% 3) Cardiac Binning: assign each readout to a cardiac phase
binOffset = 0;  % phase offset (in readouts) applied before binning
disp(['Binning k-space into ', num2str(opt.nCPhases), ' cardiac phases...']);
% Outputs of binCardiacPhases:
%   cbin_est    – cardiac bin index for each of total_samples readouts
%   meanRR      – mean RR interval (s)
%   stdRR       – standard deviation of RR intervals (s)
%   triggerTime – times of detected triggers (s)
%   triggers    – sample indices of detected triggers
%   card_int    – interpolated cardiac signal at trigger times
[cbin_est, meanRR, stdRR, triggerTime, triggers, card_int] = ...
    binCardiacPhases( ...
        total_samples, card, ...      % number of readouts, cardiac trace
        opt.TR_SG, opt.nCPhases, ...   % SG repetition time, # of bins
        opt.sgI, binOffset);           % SG interval, phase offset

% Interpolate trigger amplitudes for plotting
trigger_pks = interp1(SG_timeFull, card, triggers);

% Format and save the cardiac trigger plot
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0.5 1 0.5]);
fig_card = gcf;
saveas(fig_card, fullfile(folder_out_path, 'cardiac_triggers.png'));



% 4) Respiratory Binning: assign each readout to a respiratory phase
fix_eff = 1;  % use fixed-efficiency binning if 1
% Outputs of binRespHard:
%   rbin_est  – respiratory bin index per readout
%   resp      – (optionally updated) respiratory trace
%   bin_lines – bin threshold levels for plotting
[rbin_est, resp, bin_lines] = ...
    binRespHard( ...
        total_samples, resp, opt.sgI, ...  % # readouts, resp trace, SG interval
        opt.nRPhases, SG_timeFull, ...     % # resp bins, SG time vector
        timeFull, fix_eff);                % full time vector, fixed-efficiency flag

% Format and save the respiratory distribution plot
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0.05 1 0.5]);
fig_resp = gcf;
saveas(fig_resp, fullfile(folder_out_path, 'resp_dist.png'));



% 5) Interpolate SG-based signals onto full acquisition timeline
%   Allows plotting and analysis at every readout timepoint
resp_interp = interp1(SG_timeFull, resp, timeFull, 'spline');
card_interp = interp1(SG_timeFull, card, timeFull, 'spline');

% 6) Remove self‐gating readouts to avoid bias from oversampled SG samples
PEInd = PEInd_SG;                % full PE index list
PEInd(SG_locs_ind, :) = [];      % drop SG‐based lines
total_samples = length(PEInd);

% Remove corresponding entries in all data arrays
kspc_ds(:, :, SG_locs_ind)     = [];
cbin_est(SG_locs_ind)          = [];
rbin_est(SG_locs_ind)          = [];
resp_interp(SG_locs_ind)       = [];
card_interp(SG_locs_ind)       = [];
timeFull(SG_locs_ind)          = [];



% 7) Initialize Weighting Vectors & Compute Acceleration Rate -------

[w_init, total_unique_pairs, accel] = initWeighting( ...
PEInd, cbin_est, rbin_est, total_samples, opt, mat_size);

% 8) Crop Readouts, Coil Combine & Compute k-Space Indices -----------
[kspc_ds, PE_lin_ind, y_ind, mat_size] = cropAndIndexKspace(kspc_ds, PEInd, mat_size, opt);
dims = [mat_size(1:3),mat_size(5:6)]; % Updated Image dimensions after Cropping
n_dims = numel(dims); % Number of image dimensions
cent_slc = round(mat_size./2); % Extracting central slice for display gifs, you can select any other slice as desired


%% Recon Inits and coil combine
disp("Image recon initializations...")
y = gpuArray(single(permute(kspc_ds,[1 3 2]))); % k-space re-arranged for reconstruction
x = accumarray(y_ind,y(:),[prod(mat_size(1:4)) 1]); % Accumulating vectorized k-space into k-space locations
wo = prior_per(w_init,opt.sg_factor,opt.out_factor); % SG-based Prior
samp_scale = accumarray(y_ind,1,[prod(mat_size(1:4)) 1]); % Sampling frequency of each sampled k-space location
samp_scale(samp_scale==0) = inf; % To set unsampled locations as 0 in next line
x = x./samp_scale; % Averaging accumulated k-space according to sampling frequency
x = ifft3_shift(reshape(x,mat_size(1:4))); % Time-averaged image estimate
disp("Coil sensitivities estimation started...")
[sens,x0] = WalshCoilCombine3D(gather(x)); 
disp("Coil sensitivities estimation completed...")
opt.scale = 0.1 * max(abs(kspc_ds(:))); % Scaling factor for k-space
y = y./opt.scale; % Scaling k-space to have max value 0.1
opt.sig_scl = opt.sig./opt.scale; % Scaling noise st. dev

%%  EMORe Reconstruction
[x_em, w_em, cost_em, percentage_energy_class] = runEMOReRecon(y, x0, sens, w_init, wo, y_ind, ...
    PE_lin_ind, dims, n_dims, mat_size,cent_slc, folder_opt_path, folder_gif_path, folder_out_path, ...
    timeFull, resp_interp, card_interp, opt);

%% CS Reconstruction
if(opt.cs == 1)
[x_cs, cost_cs] = runCSRecon(y, x0, sens, w_init, y_ind, ...
    PE_lin_ind, dims, n_dims, mat_size,cent_slc, folder_opt_path, folder_gif_path, opt);
end
  


