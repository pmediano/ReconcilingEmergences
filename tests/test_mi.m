%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Unit tests for mutual information functions. Can be run from the repository's
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

%% Discrete mutual information
% Test that self-MI is equal to entropy
X = randi(2, [1000,1]);
assert(DiscreteMI(X, X) > (0.99));

% Test that MI between independent variables is close to zero
X = randi(2, [10000,1]);
Y = randi(2, [10000,1]);
assert(abs(DiscreteMI(X, Y)) < 0.001);

% Same test, but with non-binary variables
X = randi(3, [50000,1]);
Y = randi(4, [50000,1]);
assert(abs(DiscreteMI(X, Y)) < 0.01);

% Compare with analytical MI in a binary symmetric channel
p_vec = [0.1, 0.3, 0.5];
for p=p_vec
  X = randi(2, [50000,1]) - 1;
  Y = xor(X, rand(size(X)) < p);
  true_MI = 1 + p*log2(p) + (1-p)*log2(1-p);
  assert(abs(DiscreteMI(X, Y) - true_MI) < 0.01);
end


%% Gaussian mutual information
% Test that MI between independent variables is close to zero
X = randn(10000,1);
Y = randn(10000,1);
assert(abs(GaussianMI(X, Y)) < 0.001);

% Compare with analytical MI
rho_vec = [0, 0.3, 0.5, 0.8];
for rho=rho_vec
  X = mvnrnd([0,0], [1, rho; rho, 1], 50000);
  true_MI = -0.5*log(1 - rho*rho);
  assert(abs(GaussianMI(X(:,1), X(:,2)) - true_MI) < 0.01);
end

