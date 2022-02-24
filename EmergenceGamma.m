function [ gamma ] = EmergenceGamma(X, V, tau, method)
%% EMERGENCEGAMMA Compute causal decoupling criterion from data
%
%     GAMMA = EMERGENCEGAMMA(X, V) computes the causal decoupling criterion
%     gamma for the system with micro time series X and macro time series Y.
%     Micro data must be of size TxD and macro data of size TxR, where T is the
%     length of the time series and D,R the dimensions of X,V respectively.
%     Note that for the theory to hold V has to be a (possibly stochastic)
%     function of X.
%
%     GAMMA = EMERGENCEGAMMA(X, V, TAU) uses a time delay of TAU samples to
%     compute time-delayed mutual information. (default: 1)
%
%     GAMMA = EMERGENCEGAMMA(X, V, TAU, METHOD) estimates mutual info using a
%     particular METHOD. Can be 'discrete' or 'gaussian'. If empty, it will
%     try to infer the most suitable method.
%
% Reference:
%     Rosas*, Mediano*, et al. (2020). Reconciling emergences: An
%     information-theoretic approach to identify causal emergence in
%     multivariate data. https://arxiv.org/abs/2004.08220
%
% Pedro Mediano and Fernando Rosas, Aug 2020

%% Parameter checks and initialisation
if ~ismatrix(V) || ~ismatrix(X)
  error("X and V have to be 2D matrices.");
end
if size(V,1) ~= size(X,1)
  error("X and V must have the same height.");
end
if nargin < 3 || isempty(tau)
  tau = 1;
end
if nargin < 4 || isempty(method)
  if exist('OCTAVE_VERSION', 'builtin')
    isdiscrete =  (sum(abs(X(:) - round(X(:)))) + sum(abs(V(:) - round(V(:))))) < 1e-10;
  else
    isdiscrete =  iscategorical(X) || (sum(abs(X(:) - round(X(:)))) + sum(abs(V(:) - round(V(:))))) < 1e-10;
  end
  if isdiscrete
    method = 'discrete';
  else
    method = 'gaussian';
  end
end

if strcmp(lower(method), 'gaussian')
  MI_fun = @GaussianMI;
elseif strcmp(lower(method), 'discrete')
  MI_fun = @DiscreteMI;
else
  error("Unknown method. Implemented options are 'gaussian' and 'discrete'.");
end


%% Compute mutual infos and gamma
gamma = max(arrayfun(@(j) MI_fun(V(1:end-tau,:), X(1+tau:end,j)), 1:size(X,2)));

