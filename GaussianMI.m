function [ mi ] = GaussianMI(X, Y)
%%GAUSSIANMI Mutual information between Gaussian variables
%     I = GAUSSIANMI(X, Y) calculates the mutual information between vectors X
%     and Y, assuming they are jointly Gaussian distributed.
%
% Pedro Mediano and Fernando Rosas, Aug 2020

%% Parameter checks
if ~(isvector(X) && isvector(Y) && length(X) == length(Y))
  error("X and Y must be vectors of the same length.");
end

%% Compute correlation and MI
rho = corr(X(:), Y(:));
mi = -0.5*log(1 - rho*rho);

