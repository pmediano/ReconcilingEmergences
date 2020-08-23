%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Unit tests for causal emergence functions. Can be run from the repository's
% root folder with runtests('tests/').
%
% Note: these tests are stochastic, so if any of them fails try running them
% again :]
%
% Pedro Mediano and Fernando Rosas, Aug 2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add root folder to path
addpath ..

%% Test with independent Gaussian variables
X = randn([50000,2]);
V = randn([50000,1]);
assert(EmergencePsi(X, V)   < 0.01);
assert(EmergenceDelta(X, V) < 0.01);
assert(EmergenceGamma(X, V) < 0.01);

% Same, but with a different time-scale
assert(EmergencePsi(X, V, 5)   < 0.01);
assert(EmergenceDelta(X, V, 5) < 0.01);
assert(EmergenceGamma(X, V, 5) < 0.01);

% Test micro and macro MI too
[psi, v_mi, x_mi] = EmergencePsi(X, V);
assert(abs(psi) < 0.01);
assert(abs(v_mi) < 0.01);
assert(abs(x_mi) < 0.01);


%% Test with discrete downward causation example
T = 50000;
X = zeros([T,2]);
for t=2:T
  X(t,1) = xor(X(t-1,1), X(t-1,2));
  X(t,2) = rand < 0.5;
end
V = xor(X(:,1), X(:,2));

[psi, v_mi, x_mi] = EmergencePsi(X, V);
assert(abs(psi) < 0.01);
assert(abs(v_mi) < 0.01);
assert(abs(x_mi) < 0.01);

[delta, v_mi, x_mi] = EmergenceDelta(X, V);
assert(abs(delta - 1) < 0.01);
assert(all(abs(v_mi - x_mi - [1,0]) < 0.01));

assert(abs(EmergenceGamma(X, V) - 1) < 0.01);


%% Test with discrete causal decoupling example
T = 50000;
gamma = 0.99;
M = 1 + gamma*log2(gamma) + (1-gamma)*log2(1-gamma);
X = zeros([T,2]);
for t=2:T
  p = xor(xor(X(t-1,1), X(t-1,2)), rand > gamma);
  X(t,1) = rand < 0.5;
  X(t,2) = xor(p, X(t,1));
end
V = xor(X(:,1), X(:,2));

[psi, v_mi, x_mi] = EmergencePsi(X, V);
assert(abs(psi - M) < 0.01);
assert(abs(v_mi - M) < 0.01);
assert(abs(x_mi) < 0.01);

[delta, v_mi, x_mi] = EmergenceDelta(X, V);
assert(abs(delta) < 0.01);
assert(all(abs(v_mi - x_mi) < 0.01));

assert(abs(EmergenceGamma(X, V)) < 0.01);

