function [pval, permMAPDiff] = randomtest(X,Y)
% RANDOMTEST returns the p-value of the randomization test on MAP values
% 	pval = randomtest(X,Y) returns the p-value of the randomization test
% 	with the null hypothesis that the average precision values in vectors
% 	X and Y follows the same distribution. Optionally it outputs the
% 	randomization values.
%
% Copyright Sohan Seth sohan.seth@hiit.fi

origMAPDiff = mean(X) - mean(Y);
if length(X) ~= length(Y)
	error('X and Y must be of same length')
end
nSam = length(X);
nPerm = 10000; % number of permutation
X = X(:); Y = Y(:);
permMAPDiff = zeros(nPerm, 1);
for countPerm = 1:nPerm
	perm = rand(nSam, 1) > 0.5;
	permMAPDiff(countPerm) = mean([X(perm); Y(~perm)]) - mean([X(~perm); Y(perm)]);
end
pval = mean(abs(origMAPDiff) < abs(permMAPDiff));