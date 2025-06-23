function [entropyVals, bestIdx] = selectComponentByEntropy(PCIC, TR_SG)
%SELECTCOMPONENTBYENTROPY  Choose the component with lowest RR-interval entropy
%
%   [entropyVals, bestIdx] = selectComponentByEntropy(PCIC, TR_SG)
%   • PCIC:    [T × C] matrix of concatenated PCA/ICA time-series
%   • TR_SG:   SG sampling interval in seconds
%
%   Outputs:
%     entropyVals – C×1 vector of Shannon entropy values of RR intervals
%     bestIdx     – index (1..C) of component with minimum entropy
%
%   The function:
%     1) For each component, detects peaks in -gradient(signal)
%     2) Converts peak locations to RR intervals (Δt)
%     3) Filters out physiologically invalid RR (<0.35 s)
%     4) Computes histogram-based PDF over edges [0.25,2] s
%     5) Calculates entropy of the PDF
%     6) Returns the component index with lowest entropy (best rhythmicity)

% Define RR histogram edges
edges = linspace(0.25, 2, 11);
C     = size(PCIC,2);
entropyVals = zeros(C,1);

for i = 1:C
    % Compute positive gradient peaks
    sig = PCIC(:,i);
    tmp = -gradient(sig);
    tmp(tmp < 0) = 0;
    cutoff = prctile(tmp,95) * 0.5;
    tmp(tmp < cutoff) = 0;
    [~, loc] = findpeaks(tmp);

    % Convert to RR intervals
    tPks = loc * TR_SG;
    rr   = diff(tPks);
    rr(rr < 0.35) = [];

    if isempty(rr)
        entropyVals(i) = Inf;
        continue;
    end

    % PDF and Shannon entropy
    counts  = histcounts(rr, edges);
    pdfVals = counts / sum(counts);
    entropyVals(i) = -sum(pdfVals(pdfVals>0) .* log2(pdfVals(pdfVals>0)));
end

% Select component with minimum entropy
di = find(entropyVals == min(entropyVals), 1, 'first');
bestIdx = di;
end
