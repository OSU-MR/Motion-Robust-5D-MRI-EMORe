function plotEntropyHistograms(PCIC, TR_SG, folder_opt_path)
%PLOTENTROPYHISTOGRAMS  Plot RR-interval histograms and entropy/std metrics
%
%   plotEntropyHistograms(PCIC, TR_SG, folder_opt_path)
%   • PCIC:      [T x C] matrix of concatenated PCA/ICA signals
%   • TR_SG:     sampling interval of SG in seconds
%   • folder_opt_path: directory to save the histogram figure
%
%   Computes RR intervals for each component, derives entropy and std,
%   then plots a tiled histogram figure and saves it.

% Define histogram edges for RR intervals (s)
edges = linspace(0.25, 2, 11);
numComp = size(PCIC, 2);
entropyVals = zeros(numComp,1);
stdRR      = zeros(numComp,1);

% Create invisible figure for batch saving
fig = figure('Visible', 'off');

% Determine subplot layout
cols = 5;
rows = ceil(numComp/cols);

for i = 1:numComp
    % Compute gradient-based peaks
    tmp = -gradient(PCIC(:,i));
    tmp(tmp < 0) = 0;
    cutoff = prctile(tmp,95) * 0.5;
    tmp(tmp < cutoff) = 0;
    [~, loc] = findpeaks(tmp);

    % Convert peak indices to times and compute RR intervals
    timePks = loc * TR_SG;
    rr = diff(timePks);
    rr(rr < 0.35) = [];

    % Compute histogram-based PDF and entropy
    histRR = histcounts(rr, edges);
    pdfRR  = histRR / sum(histRR);
    entropyVals(i) = -sum(pdfRR(pdfRR>0) .* log2(pdfRR(pdfRR>0)));
    stdRR(i)      = std(rr);

    % Plot histogram
    subplot(rows, cols, i);
    histogram(rr, edges, 'EdgeColor', 'none');
    title(sprintf('E=%.2f S=%.2f', entropyVals(i), stdRR(i)));
    xlabel('RR interval (s)');
    ylabel('Count');
end

grid on;

% Save the figure
saveas(fig, fullfile(folder_opt_path, 'entropy_histograms.png'));
close(fig);
end
