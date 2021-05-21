function [V, ind, k, r] = newton(ker, X, tol, n)
% NEWTON Computes the Newton basis with a greedy Power Function max.
%   [V, ind, k, r] = NEWTON(ker, X, tol, n) constructs the value matrix V 
%   on X of the first n Newton basis functions. The algorithm selects k  
%   points x in X by a greedy max. of the L_infty norm of the Power 
%   on X. The algorithm stops if tol is reached before selecting n bases.
%   The values of the Power Function are in r.
%
%   NOTE: this function is a modified version of the original function by
%   Stefan Mueller and Robert Schaback, see
%       http://num.math.uni-goettingen.de/schaback/research/papers/ANBfKS_MATLAB.tgz
%

m = size(X, 1);
V = zeros(m, n);
ind = zeros(n, 1);

z = repmat(ker(0, 0), length(X), 1);
[normu, ind(1)] = max(z);
x = X(ind(1), :);
V(:, 1) = ker(X, x) / sqrt(normu);
w = V(:, 1) .^ 2;
r = zeros(n, 1);
r(1) = 1;

for k = 2 : n
    
    [maxPf, ind(k)] = max(z - w);
    r(k) = maxPf;
    if maxPf < tol
        break
    end
    x = X(ind(k), :);
    u = ker(X, x) - V(:, 1 : k - 1) * V(ind(k), 1 : k - 1)';
    normu = sqrt(z(ind(k)) - w(ind(k)));
    V(:, k) = u / normu;
    w = w + V(:, k) .^ 2;
end
ind = ind(1 : k);
ind0 = setdiff(1 : m, ind);
ind = [ind; ind0'];
V = V(ind, 1 : k);
r = sqrt(r(1 : k));

end
