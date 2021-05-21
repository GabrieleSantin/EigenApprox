function [P, n1, xx, yy, ind, con] = getPoints(omegaType, n, pointType, plotFlag)
% GETPOINTS Constructs a set of points in R^2.
%   [P, n1, xx, yy, ind] = GETPOINTS(domainType, n, pointType, plotFlag)
%   returns a set of n1 >= n points of type pointType in the domain
%   domainType and plots them in figure(plotFlag) if plotFlag ~= 0.
%   P is a n1x2 vector of points, corresponding to the points with indexes
%   ind in the grid [xx, yy].

nGuess = floor(sqrt(n));

% Generate the domain
switch omegaType
    case 'square'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        con  =  @(p) true(size(p,1),1);
    case 'disk'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        nGuess = nGuess*sqrt(4/pi);
        con  =  @(p) p(:,1).^2 + p(:,2).^2 <=  1;
    case 'pacman'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        nGuess = nGuess*sqrt(16/3/pi);
        c =  @(p) p(:,1)./sqrt(p(:,1).^2 + p(:,2).^2);
        s =  @(p) p(:,2)./sqrt(p(:,1).^2 + p(:,2).^2);
        con  =  @(p) (p(:,1).^2 + p(:,2).^2 <=  1) & (sign(c(p))+sign(s(p)) >= 0);
    case 'lune'
        xmin = 0; xmax = 1; ymin = 0; ymax = 1;
        nGuess = nGuess*sqrt(8/(pi+2));
        con  =  @(p) ((p(:,1)-0.5).^2 + (p(:,2)-0.5).^2 <=  0.25) & ...
            (p(:,1).^2 + p(:,2).^2 >=  0.25);
    case 'dbubble'
        xmin = -2; xmax = 2; ymin = -2; ymax = 2;
        nGuess = nGuess;
        con  =  @(p) ((p(:,1)-sqrt(2)/2).^2 + p(:,2).^2 <=  1) | ...
            ((p(:,1)+sqrt(2)/2).^2 + p(:,2).^2 <=  1);
    case 'lens'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        nGuess = nGuess;
        con  =  @(p) ((p(:,1)+sqrt(2)/2).^2 + p(:,2).^2 <=  1) & ...
            ((p(:,1)-sqrt(2)/2).^2 + p(:,2).^2 <=  1);
    case 'cardioid'
        xmin = -1.5; xmax = 1.5; ymin = -1.5; ymax = 1.5;
        nGuess = nGuess;
        con  =  @(p) ((p(:,1)-1).^2+p(:,2).^2+(p(:,1)-1)).^2-((p(:,1)-1).^2+p(:,2).^2) <= 0;
    case 'triangle'
        xmin = 0; xmax = 1; ymin = 0; ymax = 1;
        nGuess = nGuess*2;
        con  =  @(p) (p(:,2) >=  0) & (-p(:,1)-p(:,2)+1 >= 0) & (p(:,1)-p(:,2) >= 0);
    case 'ring'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        r2 = 0.1;
        nGuess = nGuess*20/sqrt(19*pi);
        con  =  @(p) (p(:,1).^2 + p(:,2).^2 <=  1) & (p(:,1).^2 + p(:,2).^2 >=  r2);
    case 'doublering'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        r2 = 0.05;
        nGuess = nGuess*20/sqrt(19*pi);
        con  =  @(p) (p(:,1).^2 + p(:,2).^2 <=  1) & ...
            ((p(:,1) - 0.5).^2 + p(:,2).^2 >=  r2) & ...
            ((p(:,1) + 0.5).^2 + p(:,2).^2 >=  r2);
    case 'holes'
        xmin = -1; xmax = 1; ymin = -1; ymax = 1;
        r2 = 0.02;
        nH = 16;
        c = 2 * rand(nH, 2) - ones(nH, 2);
        str = '@(p, c, r2) 1';
        for j = 1 : size(c, 1)
            str = [str, ' & ((p(:,1) - ' num2str(c(j, 1)) ...
                ').^2 + (p(:,2) - ' num2str(c(j, 2)) ').^2 >= '...
                num2str(r2) ')'];
        end
        con = str2func(str);
    case 'strip'
        xmin = -1; xmax = 1; ymin = -0.05; ymax = 0.05;
        con  = @(p) true(size(p,1),1);
    otherwise
        error(['Domain "' omegaType '" not avaiable'])
end

%Select the points inside the domain
n1 = 0;
gridSize = nGuess;

while n1<n
    switch pointType
        case 'u'
            [xx,yy] = ndgrid(linspace(xmin,xmax,gridSize),linspace(ymin,ymax,gridSize));
        case 'r'
            xx  =  rand(floor(gridSize))*(xmax-xmin)+xmin;
            yy  =  rand(floor(gridSize))*(ymax-ymin)+ymin;
        otherwise
            error(['Point type "' pointType '" not avaiable'])
    end
    grid = [xx(:),yy(:)];
    ind  =  con(grid);
    n1 = sum(double(ind));
    gridSize = gridSize+1;
end

P = [xx(ind) yy(ind)];

%Plotting
if plotFlag
    figure(plotFlag), plot(P(:,1), P(:,2), 'g.'), axis equal,
    % title(['n  =  ' num2str(size(P,1))])
end
