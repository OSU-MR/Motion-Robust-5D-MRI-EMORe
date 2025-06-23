%% Options for initialization

function [opt] = init()

% Filename and scanner version
% ************************************************
opt.fName   = [];
opt.version = 'VE';

% Read raw data from scanner (readWrapper options)
% ************************************************
opt.Nc      = [100, 0.5];   % Number of compressed coils
opt.fd      = [0.0, 0.0];   % Fraction of FE discarded
opt.nAmp    = 1;            % Noise amplification factor
opt.nStd    = 1e-5;         % Target noise std
opt.nNrm    = 0;            % Use the same noise std for each coil? 1:yes, 0:no
opt.dis     = 1;            % Display time-averaged images?
opt.samp    = [];           % Sampling pattern
opt.sDim    = [];           % Dimensions along which subsampling takes place
opt.sgAvg   = 1;            % To average the 'seg' dimension or not
opt.mxFlg   = 0;            % To find Mawell correction maps or not
opt.kNoise  = 0;            % Ignore prescan and use k-space for noise power estimation
opt.tap     = 0;            % Taper the boundaries along readout
opt.peSft   = 0;            % Circular shifting in PE direction

% Extracting self-gating signals
% **************************************************
opt.sgI     = 9;            % Self-gating interval: 0 means no self-gating
opt.sigExt  = 'PCA';        % Method to extract resp/cardiac signals. 'PCA', 'SSAFARI', or 'PT'
opt.ECG     = 0;            % Use ECG for cardiac signal?
opt.nCoils  = 12;           % # coil elements used for PCA (not used right now!)
opt.rampSS  = 25;           % Trim # self-gating samples during steady state ramp
opt.rLPF    = [0.5];        % FIR low pass filter to select respiratory signal (Hz)
opt.cBPF    = [0.5 3];      % FIR band pass filter to select cardiac signal (Hz)
opt.fOrd    = 1001;         % FIR filter order (length+1)

% Binning and data array construction
% **************************************************
opt.nRPhases    = 1;        % Number of respiratory phases to bin
opt.viewSharing = 0;        % View sharing to increase # of cardiac phases
opt.rSoft       = 1;        % Soft-gating (weighting) (1) or not (0)
opt.expOrd      = 4;        % Order of exponential for soft-gate weighting function e.g. higher = sharper cutoff
opt.nRcnR       = 1;        % Number of respiratory phases to reconstruct
opt.Rate_Max    = 23;       % Maximum allowed acceleration rate (not including repeat samples)
opt.nCPhases    = 20;       % Number of cardiac phases to bin
opt.pRcnC       = [1, 20];  % Cardiac phases to reconstruct
opt.arrhyRej    = 3;        % Reject heartbeats this many stdevs away from mean RR interval
opt.binOffset   = 0;        % Circularly shift cardiac bins
opt.recTime     = 300; %300     % Only reconstruct "recTime" worth of data in seconds, inf = all
opt.recStart    = 'center'; % 'start', 'center', 'end': where to begin trimming for "recTime"
opt.respEff     = 50;

% Misc
% ****************************************************
opt.disp = 1;
% ReVEAL4D Reconstruction
% **************************************************

% Cropping 


end

