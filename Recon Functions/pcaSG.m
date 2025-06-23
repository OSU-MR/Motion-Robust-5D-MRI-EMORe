function [PC, S] = pcaSG(A,k)
% Principle component analysis
% INPUTS:   A ~ [Time x Concatenated RO's]
%           k ~ # of components
% OUTPUTS:  PC ~ k singular vectors (principle components) with largest 
%                singular value, weighted by singular value

% Eigen-decomposition of covariance matrix
[V,D]   = eig(A'*A);

% Sort eigenvalues and eigenvectors in descending order
[~,ind] = sort(diag(D),'descend');

S = diag(D); S = S(ind);

% Sorted eigenvectors
Vs = V(:,ind);

% Singular vectors
US  = A*Vs;

% Choose first 'k' principle components
PC = US(:,1:k); 
end

