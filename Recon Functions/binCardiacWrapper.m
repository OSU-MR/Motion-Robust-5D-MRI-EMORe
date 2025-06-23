function [cBins, meanRR, stdRR, triggerTime, triggers] = binCardiacWrapper(lengthDat, C, opt)

method = opt.sigExt;


switch method
    % If PCA was used...
    case 'PCA'
        if opt.ECG
            [cBins, meanRR, stdRR, triggerTime] = binCardiacPhases_ECG(lengthDat, C, opt);            
        else
            [cBins, meanRR, stdRR, triggerTime,triggers] = binCardiacPhases(lengthDat, C, opt);
        end
        
    case 'PT'
        if opt.ECG
            [cBins, meanRR, stdRR, triggerTime] = binCardiacPhases_ECG(lengthDat, C, opt);            
        else
            [cBins, meanRR, stdRR, triggerTime,triggers] = binCardiacPhases(lengthDat, C, opt);
        end
        
    % If SSA-FARI was used...    
    case 'SSAFARI'
        cBins = binCardiacPhases_SSAFARI(lengthDat, C, opt);
        meanRR = [];
        stdRR = [];
        triggerTime = [];
        
    case 'PTSola'
        [cBins, meanRR, stdRR, triggerTime,triggers] = binCardiacPhases(lengthDat, C, opt);

end

end

