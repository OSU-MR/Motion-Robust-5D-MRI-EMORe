function [component_rank] = select_cardiac_component(signals)
% Commented out by murtaza
% ----------------------------------------------------------------%
% Preprocessing to remove large outliers in each signal
for i = 1 : size(signals, 2)
    tmp = signals(:,i);
    lo = prctile(tmp, 1);
    hi = prctile(tmp, 99);
    tmp(tmp < lo); tmp(tmp > hi) = 0;
    signals(:,i) = tmp;
end
% ----------------------------------------------------------------%
% Commented out by murtaza

% Normalize  mean and dynamic range
for i = 1 : size(signals, 2)
    tmp = rescale(signals(:,i), -1, 1);
    tmp = tmp - mean(tmp);
    signals(:,i) = tmp;
end

% Test 1
%---------------------------------------------------------------------------------------------------
spectra = fft(signals, [], 1);
Cov = diag(spectra'*spectra);

[~,component_rank] = sort(abs(Cov), 'descend');


place_holder = [];


component_index = [];



end