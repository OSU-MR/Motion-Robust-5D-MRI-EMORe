function [cBins, meanRR, stdRR, triggerTime, triggers,C] = binCardiacPhases_nofig(lengthDat,C,TR,arrhyRej,nCPhases,SG_freq,binOffset)



% C = -C;
% C = orientCardiac_2( C, TR);
arrhyRej_fac = 3;


interp_factor = 2;
if interp_factor ~= 1
    xq = (1/interp_factor):(1/interp_factor):length(C);
    C_interp = interp1(C, xq, 'spline');
end

dC = -gradient(C_interp);
% Threshold values less than 0
dC(dC < 0) = 0;
cutOff = prctile(dC, 95)*0.5;
dC(dC < cutOff) = 0;
% Find peaks of gradient - corresponds to peak systole/ejection
[pks, loc] = findpeaks(dC);



% nPeaks = length(pks);

% [pks_2, loc_2] = findpeaks(C_interp);

dC_tmp = -gradient(C_interp);
% Trim first peak for robustness
loc = loc(2:end);
loc_temp = zeros(length(loc),1);
loc_true = zeros(length(loc),1);
for i = 1 : length(loc)
    counter = 0;
    stop = 0;
    while ~stop
        counter = counter + 1;
        if (loc(i) - counter) < 1
            stop = 1;
        else
            peak(1) = dC_tmp(loc(i) - counter + 1);
            peak(2) = dC_tmp(loc(i) - counter);
 
            if counter >= 3
                if counter == 3
                    previous_location = loc(i)-counter;
                    previous_dif = inf;
                else
                    previous_location = current_location;
                    previous_dif = current_dif;
                end
                
                % Linear regression
                m = dC_tmp(loc(i)-counter+1)-dC_tmp(loc(i)-counter);
                b = dC_tmp(loc(i)-counter) - m*(loc(i)-counter);
                loc_temp(i) = loc(i) - counter + 1;
                current_location = -b/m;
                current_dif = current_location - previous_location;
                
                if abs(current_dif) >= abs(previous_dif)
                    loc_true(i) = previous_location;
                    clear peak;
                    stop = 1;
                end
            end
        end
    end
end

% loc_temp(loc_temp == 0) = [];



% loc = loc_temp;
% nPeaks = length(loc);
% nPeaks = length(loc_true);


% Time
timeSG = [1:1:length(C_interp)]*TR*(1/interp_factor); %What is interp factor?
% timePeaks = timeSG(loc);
timePeaks_true = interp1([1:length(timeSG)], timeSG, loc_true);
timePeaks = timePeaks_true;
timePeaks(isnan(timePeaks)) = [];
nPeaks = length(timePeaks);

tmp = [0:1:length(C_interp)-1]*TR*(1/interp_factor);
tmp_true = interp1([1:length(tmp)], tmp, loc_true);
triggers = tmp_true;


% Find variability statistics for arrhythmia rejection
RR = zeros(1,nPeaks-1);
for i = 1 : nPeaks-1
    RR(i) = timePeaks(i+1) - timePeaks(i); 
end

% savefig('rrhist.fig');
% saveas(gcf,'rrhist.jpg');
% Remove outliers
tmp_RR = RR;
if(arrhyRej)
tmp_RR(abs((tmp_RR-mean(tmp_RR))/std(tmp_RR)) > arrhyRej_fac) = [];
end

meanRR = mean(tmp_RR); 
stdRR = std(tmp_RR);
% medRR = median(tmp_RR); 
% percentile_95 = prctile(tmp_RR, 95);
% percentile_25 = prctile(tmp_RR, 25);

binWidth = meanRR/nCPhases;
% binWidth = percentile_25/opt.nCPhases;
% nPhases_systole = nCPhases;
nPhases_systole = ceil(nCPhases/3);

nPhases_diastole = nCPhases - nPhases_systole;


% Define which are systolic bins:

binsS = zeros(1,nCPhases);
% binsS(1:sBefore) = 1; binsS(end-sAfter+1:end) = 1; 
binsS(1:nPhases_systole) = 1;


bin = zeros(nCPhases.*(nPeaks-1),2);
% Create bins.
counter = 0;
for i = 1 : nPeaks-1
    dT = timePeaks(i+1) - timePeaks(i);
    dScaleWidth = (dT - (nPhases_systole*binWidth)) / nPhases_diastole;
    diastoleWidth(i) = dScaleWidth;
    for n = 1 : nCPhases
        counter = counter + 1;
        % Arrhythmia rejection (simple Z-score threshold)
        if ((abs((dT-meanRR) / stdRR) > arrhyRej_fac) || (dT < 0.35)) && arrhyRej
            bin(counter, 1) = NaN;
            bin(counter, 2) = NaN;
        else
            if n == 1
                binStart = timePeaks(i) - (binWidth/2) - binWidth;
                bin(counter, 1) = binStart;
                bin(counter, 2) = binStart + binWidth;
            elseif n > 1
                if binsS(n)
                    bin(counter, 1) = bin(counter-1, 2);
                    bin(counter, 2) = min(bin(counter-1, 2) + binWidth, timePeaks(i+1) - (binWidth/2) - binWidth);
                else
                    bin(counter, 1) = bin(counter-1, 2);
                    bin(counter, 2) = bin(counter-1, 2) + dScaleWidth;
                end
            end
        end
    end
end



% Now we need to divide the k-space into the bins
% Timestamp for each line
timeDat = (1:1:lengthDat)*(TR/SG_freq);


DatBin = zeros(1,length(timeDat));

for i = 1 : length(timeDat)
    for j = 1 : length(bin)
        if j == 1
            % binStart = timePeaks(1) - (binWidth/2);
%             if timeDat(i) < bin(j, 2) && timeDat(i) >= binStart
            if timeDat(i) < bin(j, 2) && timeDat(i) >= bin(j, 1)
                DatBin(i) = j;
            end
        elseif j > 1
%             if timeDat(i) < bin(j) && timeDat(i) >= bin(j-1)
            if timeDat(i) < bin(j, 2) && timeDat(i) >= bin(j, 1)
                DatBin(i) = j;
            end
        end
    end
end

% Regions to trim later
DatBin(DatBin ==0) = NaN;
% offset 1st bin into late diastole
DatBin = DatBin + binOffset;

% modulo operator by number of cardic phases
DatBin_Phases = mod(DatBin,nCPhases);
DatBin_Phases(DatBin_Phases == 0) = nCPhases;

% Output
cBins = DatBin_Phases;

triggerTime = [];

trigger_pks = interp1(timeSG, C_interp, timePeaks_true);



end

