function DM = distanceMatrix(x, y)
% DISTANCEMATRIX Constructs the distance matrix of 2 sets of points in R^d.
%   DM = DISTANCEMATRIX(x, y) constructs the distance matrix of the set
%   of points x and y, i.e.,
%       DM(i,j) = || x(i, :) - y(j, :) ||_2
%
%   NOTE: this function is a modified version of the original function by
%   Gregory Fasshauer, see
%           http://www.math.iit.edu/~fass/590/handouts/DistanceMatrix.m
%

[nx, dx] = size(x);
[ny, dy] = size(y);
if dx~=dy
    error('x and y contain points in different space dimension')
end
DM = repmat(sum(x.*x, 2), 1, ny) - 2*x*y' + repmat((sum(y.*y, 2))', nx, 1);
DM = sqrt(DM);
end