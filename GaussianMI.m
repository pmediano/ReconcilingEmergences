function [ mi ] = GaussianMI(X, Y)
%%GAUSSIANMI Mutual information between Gaussian variables
%     I = GAUSSIANMI(X, Y) calculates the mutual information between matrices X
%     and Y, assuming they are jointly Gaussian distributed. Matrices must be
%     of size [T, Dx] and [T, Dy] respectively, where Dx,Dy are the dimensions
%     of X,Y and T is the number of samples.
%
% Pedro Mediano and Fernando Rosas, Aug 2020

%% Parameter checks
if ~(ismatrix(X) && ismatrix(Y) && size(X,1) == size(Y,1))
  error("X and Y must be matrices of the same height.");
end
Dx = size(X, 2);
Dy = size(Y, 2);

%% Compute correlation and MI
rho = corr([X, Y]);
mi = 0.5*(logdet(rho(1:Dx,1:Dx)) + logdet(rho(Dx+1:end,Dx+1:end)) - logdet(rho));

end


function LD = logdet(V)
%%LOGDET Numerically safe computation of the log-determinant of posdef matrix V
%
%     This function is taken from Lionel Barnett's MVGC toolbox:
%
%     https://www.github.com/SacklerCentre/MVGC2

[L,cholp] = chol(V);
if cholp == 0
	LD = 2*sum(log(diag(L)));
else
	DV = det(V);
	if abs(imag(DV)) > sqrt(eps)
		LD = NaN;
	else % give it the benefit of the doubt...
		LD = log(real(DV));
	end
end
end

