%%
%   Approximation of the eigenbasis
%


%%
% close all
clear all
addpath('Functions/')


%% Parameters
kerType = 'gauss';  % kernel (see getRbf.m)
omega = 'disk';     % domain Omega (see getPoints.m)
mO = pi;            % Leb. measure of Omega
ep = 1;             % shape parameter
m = 10000;          % grid size
n = 50;             % subspace size
tol = 1e-10;        % tolerance for the greedy alg.
nPlot = 9;         % number of basis to show

%% Loading
ker = getRbf(kerType); % radial basis
ker = @(x, y) ker(ep, distanceMatrix(x, y)); % symmetric kernel
[X, m, xx, yy, indPlot, con] = getPoints(omega, m, 'u', 0); % starting grid


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
% set(f, 'Position', [20, 20, 550, 650]);

for j = 1 : nPlot
    subplot(3, 3, j)
    titleString = '';
    fplot = 0 * xx;
    fplot(indPlot) = Vu(:, j);
    set(gca, 'FontSize', 18)
    % plotSurf(xx, yy, indPlot, con, Vu(:, j), titleString, 1)
    contour(xx, yy, fplot, 15, 'linewidth', 2), %hold on
    axis equal
    set(gca, 'XTickLabel', '', 'YTickLabel', '')
end


% print(f, '-depsc', ['Figures/Cover_' omega '_surf.eps'])
% print(f, '-depsc', ['Figures/Cover_' omega '_contour.eps'])
% print(f, ['Figures/Cover_' omega '_surf.png'])
print(f, ['Figures/Cover_' omega '_contour.png'])

