%%
%   Decay rate of the point-based eigenvalues of the Matern kernel 
%   vs.
%   Expected approximation order in the corresponding Sobolev space
%
%   and
%
%   Approximation rate of the point-based eigenvalues of the Matern kernel 
%   vs.
%   Expected approximation rate
%
%


%%
close all 
clear all
addpath('Functions/')
printFlag = 0; % print figures?


%% Parameters
kerType = 'mat1';   % mat0, mat1, mat2, mat3
beta = 1;           % smoothness of the kernel
omega = 'disk';     % domain Omega (see getPoints.m)
mO = pi;            % Leb. measure of Omega 
d = 2;              % space dimension
ep = 1;             % shape parameter
N = 10000;          % grid size
n = 100;            % subspace size
tol = 1e-17;        % tolerance for the greedy alg.


%% Loading  
e = -2 * (beta + d) / d; % exponent of Sobolev width asympt.
ker = getRbf(kerType); % radial basis
ker = @(x, y) ker(ep, distanceMatrix(x, y)); % symmetric kernel
X = getPoints(omega, N, 'u', 0); % starting grid
N = length(X); % update N 


%% Newton basis (with greedy L_inf maximization of the Power Function)
[V, ind, n] = newton(ker, X, tol, n); 


%% Discrete eigenvalues
G = V' * V * (mO / N); % L_2 gramian matrix
l = svd(G); % eigenvalues 
N = (1 : n)' .^ e; % theoretical decay
c = mean(l ./ N); % rescaling constant


%% Plots
f1 = figure(1); 
plot(X(ind(1 : n), 1), X(ind(1 : n), 2), 'r*'), axis equal
title('Points selected by the greedy algorithm')

jj = 1 : 20 : n;
NN = N(jj);
f2 = figure(2);
set(gca, 'FontSize', 18)
semilogy(l, 'b-', 'linewidth', 2), grid on, hold on
semilogy(c * N, 'g-.', 'linewidth', 2), 
semilogy(jj, c * NN, 'go', 'linewidth', 2), hold off,
legend('Discrete eigenvalues', ...
    ['Sobolev n-width asympt. (c = ' num2str(c, '%2.2f') ' )']),
xlabel('n')


%% Difference of the true and discrete eigenvalues
discreteSum = zeros(n, 1);
for j = 1 : n
    ll = svd(G(1 : j, 1 : j));
    discreteSum(j) = sum(ll);
end


%% Estimate of the coefficient
e1 = - 2 * (beta/ 2) / d; % proven rate
e2 = - 2 * (beta + d / 2) / d; % guessed rate
N1 = (1 : n)' .^ e1; % decay 1
N2 = (1 : n)' .^ e2; %decay 2
c1 = mean((ker(0,0) * mO - discreteSum) ./ N1); % rescaling constant 1
c2 = mean((ker(0,0) * mO - discreteSum) ./ N2); % rescaling constant 2


%% Plots
NN1 = N1(jj);
NN2 = N2(jj);
f3 = figure(3);
set(gca,'FontSize', 18)
semilogy(ker(0,0) * mO - discreteSum, 'b-','linewidth', 2), grid on, hold on
semilogy(c1 * N1, 'g-.', 'linewidth', 2), 
semilogy(c2 * N2, 'r-.', 'linewidth', 2), 
semilogy(jj, c1 * NN1, 'go', 'linewidth', 2), 
semilogy(jj, c2 * NN2, 'r^', 'linewidth', 2), hold off,
legend('Eigenvalues difference', ...
     ['Theoretical approx. rate (c = ' num2str(c, '%2.2f') ' )'],...
     ['Guessed approx. rate (c2 = ' num2str(c2, '%2.2f') ' )']),
xlabel('n')


%% Print 
if printFlag 
    print(f1, '-depsc', ['Figures/points_' num2str(beta)])
    print(f2, '-depsc', ['Figures/order_' num2str(beta)])
    print(f3, '-depsc', ['Figures/difference_' num2str(beta)])
end
