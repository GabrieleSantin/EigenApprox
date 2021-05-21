%%
%   Approximation of the eigenbasis in 1D
%


%% 
close all
clear all
addpath('Functions/')


%% Parameters
kerType = 'mat0';  % kernel (see getRbf.m)
a = 1;
b = 2;
mO = b - a;         % Leb. measure of Omega
ep = 1;             % shape parameter
m = 10000;          % grid size
n = 50;             % subspace size
tol = 1e-10;        % tolerance for the greedy alg.
nPlot = 12;         % number of basis to show


%% Loading
ker = getRbf(kerType); % radial basis
ker = @(x, y) ker(ep, distanceMatrix(x, y)); % symmetric kernel
X = linspace(a, b, m)';


%% Newton basis (with std greedy alg. by L_infty maximization)
[V, ind, n] = newton(ker, X, tol, n);


%% Approximation of the eigenbasis
G = V'*V*(mO/m); % L_2 gramian matrix
[Q, L] = svd(G); % eigenbasis
l = diag(L); % eigenvalues
Vu = V * Q; % eigenbasis (normalized in the native space) evaluated on X
[~, invInd] = sort(ind, 'ascend');
Vu = Vu(invInd, :); % sort the basis according to the greedy selection

%% Plots
f = figure(1);
nPlot = min(nPlot, n);

for j = 1 : nPlot
    subplot(4, 3, j)
    titleString = ['Eigenbasis ' num2str(j)];
    fplot = Vu(:, j);
    set(gca, 'FontSize', 18)
    plot(X, fplot, 'linewidth', 2), %hold on
    set(gca, 'XTickLabel', '', 'YTickLabel', '')
end
% print(f, '-depsc', ['Copertina' omega '.eps'])


