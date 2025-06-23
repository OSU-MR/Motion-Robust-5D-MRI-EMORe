function [rBins, rWeights, respEfficiency] = binRespWrapper(lengthDat, R, opt, param)

% Check which method was used to extract signal
switch opt.sigExt
    % If PCA was used...
    case 'PCA'
    
        signal = R;
        
    case 'PT'
        
        signal = R;
        
    % If SSA-FARI was used...
    case 'SSAFARI'
        
        phase = atan2(R(:,2),R(:,1));
        signal = phase;
        
     % If Pilot Tone was used (Sola)
    case 'PTSola'
        signal = R;
        
end



switch opt.rSoft
    % Hard binning
    case 0
        switch opt.rOverlap
            case 0
                rBins(1,1,:,:) = binRespHard(lengthDat, signal, opt);
%                 rBins = permute(rBins, [1,2,4,3]);
%                 rWeights = ones(1,1,lengthDat,opt.nRPhases);
                rWeights = logical(rBins);
                respEfficiency = 1;
            case 1
                rBins(1,1,:,:) = binRespHard_Overlap(lengthDat, signal, opt);
                rWeights = logical(rBins);
                respEfficiency = 1;
        end

        
    % Soft-gating
    case 1
%         [rWeights(1,1,:), respEfficiency] = binRespSoft(lengthDat, signal, opt, param);
        [rWeights(1,1,:,:), respEfficiency] = binRespSoft_v2(lengthDat, signal, opt, param);
        rBins = ones(1,1,lengthDat);
        
end
    
        
        

% % Retain only # respiratory phases specified for reconstruction
% nRcnR = min(opt.nRcnR, max(rBins(1,1,:)));
% 
% inds = zeros(lengthDat,1,'logical');
% rBinsTmp = rBins;
% for i = 1 : nRcnR  
%     % Next most frequent respiratory phase
%     modeR = mode(rBinsTmp);
%     
%     % Update indices
%     tmp = find(rBins == modeR);
%     inds(tmp) = 1;
%     
%     % Remove current mode
%     rBinsTmp(rBinsTmp == modeR) = [];
%         
% end
% 
% rBins(inds == 0) = NaN;

end

