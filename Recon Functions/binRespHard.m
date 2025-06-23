function [rBins,R,bin_lines] = binRespHard(lengthDat, R, SG_freq,nRPhases,SG_timeFull,timeFull,equal_data_distribution)

rAmpDat = R;
rAmpDat = interp1(SG_timeFull,rAmpDat,timeFull,'linear'); 

rBins = nan(size(rAmpDat));
% Flag to choose binning method
% If true, use equal data distribution; if false, use equal amplitude range

if equal_data_distribution
    % Calculate the 25th, 50th, and 75th percentiles
    p0 = prctile(rAmpDat, 0);
    p25 = prctile(rAmpDat, 25);
    p50 = prctile(rAmpDat, 50);
    p75 = prctile(rAmpDat, 75);
    p100 = prctile(rAmpDat, 100);


    bin_lines = [p25 p50 p75];

    % Assign bin numbers based on percentile thresholds
    rBins(rAmpDat <= p25) = 4;                       % Bin 1: Values between 0% - 25%
    rBins(rAmpDat > p25 & rAmpDat <= p50) = 3;        % Bin 2: Values between 25% - 50%
    rBins(rAmpDat > p50 & rAmpDat <= p75) = 2;        % Bin 3: Values between 50% - 75%
    rBins(rAmpDat > p75) = 1;                         % Bin 4: Values between 75% - 100%
else
    % Calculate equally spaced amplitude bin edges with outlier robust estimation
    minAmp = prctile(rAmpDat, 5);  % Use 5th percentile as lower bound to avoid outliers
    maxAmp = prctile(rAmpDat, 95); % Use 95th percentile as upper bound to avoid outliers
    binEdges = linspace(minAmp, maxAmp, 5); % Four bins, hence five edges
    bin_lines = [binEdges(2:4)];

    % Assign bin numbers based on amplitude range
    rBins(rAmpDat <= binEdges(2)) = 4;                       % Bin 1: Values between minAmp - 25% range
    rBins(rAmpDat > binEdges(2) & rAmpDat <= binEdges(3)) = 3; % Bin 2: Values between 25% - 50% range
    rBins(rAmpDat > binEdges(3) & rAmpDat <= binEdges(4)) = 2; % Bin 3: Values between 50% - 75% range
    rBins(rAmpDat > binEdges(4)) = 1;                         % Bin 4: Values between 75% - maxAmp
end

% Set rBins to NaN where rAmpDat was originally NaN
% rBins(isnan(R)) = NaN;

% Plot rAmpDat with min and max lines
figure;
plot(timeFull,rAmpDat);
title('Repiratory Signal');
xlabel('Time');
ylabel('Amplitude');
hold on;
if ~equal_data_distribution
    yline(minAmp, '--r', 'Min Bound'); % Plot min bound
    yline(maxAmp, '--b', 'Max Bound'); % Plot max bound
else
    yline(p0, '--r'); % Plot min bound
    yline(p25, '--r'); % Plot min bound
    yline(p50, '--r'); % Plot min bound
    yline(p75, '--r'); % Plot min bound
    yline(p100, '--r'); % Plot min bound

end
hold off;



end
