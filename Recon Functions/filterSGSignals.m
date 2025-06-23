function [resp, card] = filterSGSignals(SG, opt, folder_opt_path)
%FILTERSGSIGNALS   Extract respiratory & cardiac signals from SG data
%
%   [resp, card] = filterSGSignals(SG, opt, folder_opt_path)
%   Inputs:
%     SG              – [RO × coils × SG_samples] self‐gating data
%     opt.rBPF        – respiratory bandpass frequencies [Hz]
%     opt.cBPF        – cardiac    bandpass frequencies [Hz]
%     opt.FS_SG       – SG sampling rate (Hz)
%     opt.TR_SG       – SG repetition time (s)
%     folder_opt_path – directory to save diagnostic figures
%
%   Outputs:
%     resp – respiratory trace (vector)
%     card – cardiac     trace (vector)

    % Build Casorati matrix: [Time × Readouts]
    caso = reshape(SG, size(SG,1)*size(SG,2), size(SG,3));
    caso = abs(caso');
    caso = caso - repmat(mean(caso,1), [size(caso,1),1]);
    disp("SG signal extracted...");
    
    %% Cardiac & Respiratory Filtering
    fOrd = 1001;
    rFilt = fir1(fOrd-1, opt.rBPF/(opt.FS_SG/2), hann(fOrd));
    cFilt = fir1(fOrd-1, opt.cBPF/(opt.FS_SG/2), hann(fOrd));
    resp_caso = zeros(size(caso));
    card_caso = zeros(size(caso));
    for j = 1 : size(caso,2)
        resp_caso(:,j) = conv(caso(:,j), rFilt, 'same');
        card_caso(:,j) = conv(caso(:,j), cFilt, 'same');
    end
    disp("SG Cardiac and Respiratory Filtering...");

    %% Extract Respiratory Signal
    resp = pcaSG(resp_caso, 1);
    [muGMM,~,~] = gmm(resp(resp>prctile(resp,5) & resp<prctile(resp,95)));
    if muGMM < mean(resp)
        resp = -resp;
        disp("Signal inverted to have end-expiratory bin on top...");
    else
        disp("Correct orientation of resp signal.....");
    end
    disp("Respiratory signal extracted...");
    
    %% Extract Cardiac Signal via helper
    card = extractCardSignal(card_caso, opt, folder_opt_path);
end
