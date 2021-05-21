function rbf = getRbf(type)
% GETRBF Returns a strictly positive definite (in R^2) RBF kernel.
%   rbf  =  GETRBF(type) returns a function handle representing the radial
%   kernel of type 'type'.
%   The handled function is in the form
%       rbf(ep, r)
%   where ep is the shape parameter and r is the radial variable.

switch type
    case 'gauss',
        rbf = @(ep, r) exp(-(ep*r).^2);
    case 'imq',
        rbf = @(ep, r) 1./sqrt(1+(ep*r).^2);
    case 'genimq',
        rbf = @(ep, r) 1./(1+(ep*r).^2).^2;
    case 'iq',
        rbf = @(ep, r) 1./(1+(ep*r).^2);
    case 'mat0',
        rbf = @(ep, r) exp(-(ep*r));
    case 'mat1',
        %         rbf = @(ep, r) exp(-(ep*r)).*(1+(ep*r));
        rbf = @(ep, r) exp(-(sqrt(3) * ep*r)).*(1 + (sqrt(3) * ep*r));
    case 'mat2',
        %         rbf = @(ep, r) exp(-(ep*r)).*(3+3*(ep*r)+(ep*r).^2);
        rbf = @(ep, r) exp(-(sqrt(5) * ep*r)).*(1 + sqrt(5) *(ep*r) + 5/3 *(ep*r).^2);
    case 'mat3',
        rbf = @(ep, r) exp(-(ep*r)).*(15+15*(ep*r)+6.*(ep*r).^2+(ep*r).^3);
    case 'laggauss1',
        rbf = @(ep, r) exp(-(ep*r).^2).*(2-(ep*r).^2);
    case 'laggauss2',
        rbf = @(ep, r) exp(-(ep*r).^2).*(3-3*(ep*r).^2+0.5*(ep*r).^4);
    case 'lingenimq'
        rbf = @(ep, r) (2-(ep*r).^2)./(1+(ep*r).^2).^4;
    case 'wen20',
        rbf = @(ep, r) max(1-(ep*r),0).^2;
    case 'wen21',
        rbf = @(ep, r) (max(1-(ep*r),0).^4).*(3*(ep*r)+1);
    case 'wen23',
        rbf = @(ep, r) (max(1-(ep*r),0).^8).*(480*(ep*r).^3+375*(ep*r).^2+120*ep*r+15);
    otherwise
        error(['RBF "' type '" not avaiable'])
end