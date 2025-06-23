function [center, stdev, phi ] = gmm( epsilon )
%% gmm(): Gaussion Mixture Model
% =========================================================================
%% Default Options
gmmComp     = 2;                                                     % Number of Gaussian components
delta       = 2;                                                     % Minimum # of std's side distributions must be away from central
p           = 0.5;                                                   % Minimum "strength" of central distribution
Kmax        = 1000;                                                 
tol         = var(epsilon)*1e-4;                                    
%% INITIALIZATION
epsilon = sort(epsilon);

% Define contraints on the EM algorithm




% Choose initial values for the parameters

% Set 'm' to the number of data points
m = size(epsilon, 1);

% Use the overall variance of the dataset as the initial variance for each cluster.
gamma = ones(1, gmmComp) * sqrt(var(epsilon))/2;

% Initialize means similar to how we expect the solution to look like
mu = [prctile(epsilon,25),prctile(epsilon,75)];

prob = [0.5 0.5];

%% EXPECTATION MAXIMIZATION

% Matrix to hold the probability that each data point belongs to each
% cluster. One row per data point, one column per cluster
W = zeros(m, gmmComp); % !!!

% Loop until convergence.
for k = 1 : Kmax
    %% Expectation
    % Calculate the probability for each data point for each distribution.

    % Matrix to hold the pdf value for every data point for every cluster.
    % One row per data point, one column per cluster.
%     pdf = zeros(m, gmmComp); %!!! uncomment later
     pdf = zeros(m, gmmComp); 
    
    % For each cluster...
    for j = 1 : gmmComp
        % Check if gamma is equal to zero since that will cause a division
        % by zero when evaluating the gaussian
        if gamma(j) <= std(epsilon)/50 
            % Set to lower limit
            gamma(j) = std(epsilon)/50;
        end
        % Evaluate the Gaussian for all data points for cluster 'j'.
        pdf(:, j) = gaussian1D(epsilon, mu(j), gamma(j));
    end       
    
    
    % Multiply each pdf value by the prior probability for each cluster.
    pdf_w = bsxfun(@times, pdf, prob);
    
    % Divide the weighted probabilities by the sum of weighted probabilities for each cluster.
    sum_pdf_w = sum(pdf_w, 2);
    
    % Make sure sum_pdf_w > 0 to prevent divide by zero error
    if min(sum_pdf_w(:)) == 0
        sum_pdf_wTMP = sum_pdf_w(:);
        secondMin = min(sum_pdf_wTMP(sum_pdf_wTMP > 0));
        sum_pdf_w(sum_pdf_w==0) = secondMin;
    end
    
    W = bsxfun(@rdivide, pdf_w, sum_pdf_w);

    %% Maximization
    % Calculate the probability for each data point for each distribution
    
    % Store the previous means so we can check for covergence.
    prevMu = mu;

    % For each of the clusters...
    for j = 1: gmmComp
        
        % Calculate the prior probability for cluster 'j'
        prob(j) = mean(W(:, j));
        
        % Calculate the new mean for cluster 'j' by taking the weighted
        % average of *all* data points.
        mu(j) = weightedAverage(W(:, j), epsilon);
        
        % Calculate the variance for cluster 'j' by taking the weighted
        % average of the squared differences from the mean for all data
        % points.
        variance = weightedAverage(W(:, j), (epsilon - mu(j)).^2);
        
        % Calculate gamma by taking the square root of the variance.
        gamma(j) = sqrt(variance);
        
    end
    %% CONSTRAINTS
    % Adjust mean, gamma, and prob if necessary given the defined contraints
    
    % Location of the left-most peak
    [~,mnInd] = min(mu);
    
    % Location of the right-most peak
    [~,mxInd] = max(mu); 
    
    % Push leftmost peak to the left if needed
    mu(mnInd) = min(mu(mnInd), mean(epsilon));
    
    % Push rightmost peak to the right if needed
    mu(mxInd) = max(mu(mxInd), mean(epsilon));
   
    
    % Check for convergence.
    
    if (abs(mu - prevMu) < tol)
        break
    end
    
end

figure;
plotGMM(epsilon,mu,gamma,prob);
xlabel('Respiratory Depth');
ylabel('Probability');

G1 = gaussian1D(epsilon, mu(1), gamma(1))*prob(1);
G2 = gaussian1D(epsilon, mu(2), gamma(2))*prob(2);

height = [max(G1),max(G2)];

% Outputs
ind = find(height == max(height));
center = mu(ind);
stdev = gamma(ind);
phi = prob(ind);

