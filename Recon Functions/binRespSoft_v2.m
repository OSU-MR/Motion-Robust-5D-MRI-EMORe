% **************************************************************************************************
%
% This is a function which implements a soft-gating approach to respiratory
% binning. From this, a single "phase" is constructed.
%
% v2: Now supports creation of multiple respiratory bins, each with its own
% weighting function
%
%
% **************************************************************************************************

function [rWeights, respEfficiency] = binRespSoft_v2(lengthDat, R, sgI,nRPhases,respEff)
% ==================================================================================================
% RESAMPLE SG TO MATCH MEASURED DATA
% --------------------------------------------------------------------------------------------------
rAmpDat = R;
rAmpDat = (resampleSINC(rAmpDat',sgI))';

% rAmpDat = interp1(rAmpDat,1+(1/sgI):(1/sgI):length(rAmpDat)+1,'linear'); 
tmp = zeros(lengthDat,1);
tmp(1:sgI-1) = ones(sgI-1,1)*R(1); 
tmp(sgI:end) = rAmpDat(1:end-sgI+1);
rAmpDat = tmp;
% --------------------------------------------------------------------------------------------------

% ==================================================================================================
% FIND CENTER OF RESPIRATORY BINS\PHASES, BY AMPLITUDE
% --------------------------------------------------------------------------------------------------

% Use Gaussian mixture model to estimate polarity of end-expiration
[muGMM, sigmaGMM, phiGMM] = gmm(R);
% flip respiratory signal so that end-expiration is positive
if muGMM < mean(R)
    R = -R;
    rAmpDat = -rAmpDat;
    muGMM = -muGMM;
end

% If only one phase is chosen, use Gaussian mixture model to find center:
if nRPhases == 1
    mu = muGMM;
    sigma = sigmaGMM;
    sigma_orig = sigma;

% If multiple phases are chosen, divide respiratory profile into equal
% segments and take center of each segment:
else
    mxR = prctile(R,95); mnR = prctile(R,5);
    dR = (mxR - mnR)/nRPhases;
    mu(1) = mxR - (dR/2);
    for i = 2 : nRPhases
        mu(i) = mu(i-1) - dR;
    end
    sigma = dR;
    sigma_orig = sigma;
end



% ==================================================================================================
% CALCULATE WEIGHTING FUNCTIONS
% --------------------------------------------------------------------------------------------------

% +++++ Defined for fixed respiratory efficiency per bin +++++

order = 4;
eff_target = respEff;
for phase = 1 : nRPhases

    % **************************************************************************
    % This bit needs work
    sigma0 = std(rAmpDat)^2; % sigma0 = std(rAmpDat)
    if nRPhases == 1
        delta = linspace(min(rAmpDat),max(rAmpDat), 100);
    else
        delta = 0;
    end
    sigma = linspace(0, 4*sigma0, 500);
    % **************************************************************************
    
    eff = zeros(length(delta),length(sigma));
    for i = 1 : length(delta)
        for j = 1 : length(sigma)
            w_tmp = exp(-((rAmpDat-(mu(phase)+delta(i))).^order)/(sigma(j)^2));
            eff(i,j) = (sum(w_tmp(:).^2) / length(w_tmp)) * 100;
        end
    end

    % for every shift delta, find sigma which produces efficiency closest to
    % desired efficiency
    sigma_tmp = zeros(1,length(delta));
    for i = 1 : length(delta)
        [~,ind] = min(abs(eff(i,:) - eff_target));
        sigma_tmp(i) = sigma(ind);
    end
    % figure; plot(delta, sigma_tmp);

    [~,ind] = min(sigma_tmp);
    delta = delta(ind);
    sigma = sigma_tmp(ind);
    sigmaFinal(phase) = sigma;
    FWHM(phase) = 2*nthroot(-sigma^2*log(0.5),order);

    % Recompute weights
    rWeights(:,phase) = exp(-((rAmpDat-(mu(phase)+delta)).^order)/(sigma^2));

    
    % Display Rate/Efficiency
    % asym_perc = (1); % asym_perc = (2/3); don't count asymmetric echo as part of acceleration rate
    % Rate = opt.repeatFrac*prod(opt.ArraySize(3:end))/(((sum(rWeights(:).^2)*asym_perc)));

    respEfficiency = (sum(rWeights(:,phase).^2) / length(rWeights(:,phase))) * 100;
    disp(['Respiratory efficiency is: ',num2str(respEfficiency),'%']);
end
sigma = sigmaFinal;
disp_fig = 1;
% =========================================================================
if disp_fig
    
    % Plot bins on profile
    figure;
    t = [1:1:length(R)] * 4.*sgI./1000;
    plot(t,R,'b'); hold on; %plot(t,R,'.k'); hold off;

    for phase = 1 : nRPhases
        binLo = mu(phase) - (sigma(phase)/2);
        binHi = mu(phase) + (sigma(phase)/2);
        FWHMLo = mu(phase) - (FWHM(phase)/2);
        FWHMHi = mu(phase) + (FWHM(phase)/2);
        
        plot([1,t(end)],[mu(phase), mu(phase)], '--k'); hold on;
        plot([1,t(end)],[FWHMLo,FWHMLo],'--r'); hold on;
        plot([1,t(end)],[FWHMHi,FWHMHi],'--r'); hold on;
%         plot(t,R,'b'); hold on; %plot(t,R,'.k'); hold off;
    end

    % Plot weights fxns on respiratory histogram
    figure;
    rAmpDat_Sort = sort(rAmpDat,'ascend');
    h = histogram(rAmpDat_Sort, 70); hold on;
    h.FaceAlpha = 0.2;
    factor = max(h.Values);
    
    for phase = 1 : nRPhases
        rWeights_Sort = exp(-((rAmpDat_Sort-mu(phase)).^order)/(sigma(phase)^2));
        plot(rAmpDat_Sort, rWeights_Sort*factor, 'LineWidth',4); hold on
        axis([-inf, inf, 0, factor+(factor/4)]);
    end
end
placeHolder = [];
% maybe some more plotting here for meetings/testing purposes
% =========================================================================

