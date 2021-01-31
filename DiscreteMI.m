function [ mi ] = DiscreteMI(X, Y)
%%DISCRETEMI Mutual information between discrete variables
%     I = DISCRETEMI(X, Y) calculates the mutual information between vectors X
%     and Y taking each distinct value as a different discrete symbol.
%
% Pedro Mediano and Fernando Rosas, Aug 2020

%% Parameter checks
if ~(isvector(X) && isvector(Y) && length(X) == length(Y))
  error("X and Y must be vectors of the same length.");
end


%% Compute marginal entropies
H_fun = @(p, tol) -sum(p(p > tol).*log2(p(p > tol)));

[~,~,sX] = unique(X, 'rows');
tx = tabulate(sX);
Hx = H_fun(tx(:,3)/100, 1e-8);

[~, ~, sY] = unique(Y, 'rows');
ty = tabulate(sY);
Hy = H_fun(ty(:,3)/100, 1e-8);


%% Make joint distribution, and estimate MI
j = max(sY)*(sX - 1) + sY;
tj = tabulate(j);
Hxy = H_fun(tj(:,3)/100, 1e-8);

mi = Hx + Hy - Hxy;

