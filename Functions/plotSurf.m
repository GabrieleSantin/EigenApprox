function varargout = plotSurf(xx, yy, ind, con, f, titleString, figN, fileName)
% PLOTSURF Plots the graph of a function R^2 -> R.

x = xx(ind);
y = yy(ind);
t = delaunay(x, y);
px = mean(x(t), 2);
py = mean(y(t), 2);
ind = (con([px py]) > 0);
t = t(ind, : );
fig = figure(figN);
trimesh(t, x, y, f, f, ...
    'facecolor', 'interp', 'edgeColor', 'k', 'edgeAlpha', 0.1)
set(gca, 'Fontsize', 14)
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
zlabel('z', 'FontSize', 14, 'Rotation', 0);
ylim(colorbar(), [min(f) max(f)]),
title(titleString)
colormap('default')

if  exist('fileName', 'var') && ~isempty(fileName)
    colormap(gray)
    print(fig, '-deps', strcat(fileName, '.eps'))
end

if nargout
    varargout{1} = fig;
end
end

