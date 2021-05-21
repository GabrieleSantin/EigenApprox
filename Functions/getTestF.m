function f = getTestF(type)
% GETTESTF Returns a test function in R^2.
%   f  =  GETTESTF(type) returns a function handle representing the test  
%   function of type 'type'. 
%   The handled function is in the form
%       f(x, y) 

switch type
    case 'franke',
        f = @(x,y) franke(x,y);
    case 'oscill',
        f = @(x,y) cos(10*(x+y));
    case 'poly',
        f = @(x,y) (x+y-1).^9;
    case 'exp',
        ff = @(x,y,xx,yy,ep) exp(-ep^2*((x-xx).^2 + (y-yy) .^2));
        f = @(x,y) ff(x,y,0,-0.9,0.5) + ff(x,y,0,1.1,0.5) ...
            + ff(x,y,-0.4,0,0.5) + ff(x,y,0.2,0,0.6);
    case 'natgauss',
        ff = @(x,y,xx,yy,ep) exp(-ep^2*((x-xx).^2+(y-yy).^2));
        f = @(x,y) ff(x,y,0,-1.2,1) + 2*ff(x,y,-0.4,0.5,1) ...
            - 2*ff(x,y,-0.4,1.1,1) + 3*ff(x,y,1.2,1.3,1);
    case 'pacman',
        f = @(x,y) exp(abs(x-y)) - 1;
    case 'trig',
        f = @(x,y) sin(x)+cos(y);
    case 'const',
        f = @(x,y) x./x;
    case 'jumps0',
        f = @(x,y) -double(x<0 & y<0) + double(x>0 & y>0);
    case 'jumps1',
        f = @(x,y) -double(x<0).*x + double(x>0).*x...
                -double(y<0).*y + double(y>0).*y;
    case 'jumps2',
        f = @(x,y) -double(x<0).*x.^2 + double(x>0).*x.^2 ...
                -double(y<0).*y.^2 + double(y>0).*y.^2;    
    case 'jumps2a',
        f = @(x,y) -double(x<0).*x + double(x>0).*x.^2 ...
                -double(y<0).*y.^2 + double(y>0).*y.^2;   
    case 'kes',
        f = @(x,y) kes(x,y);
    otherwise
        error(['Function "' type '" not avaiable'])
end
end


function z = kes(xx,yy)
% Note: From Richard Franke's test data, see
%   https://www.staff.uni-giessen.de/odavydov/tsfit/Franke_test_data/README

[m,n] = size(xx);
N = m*n;
x = reshape(xx,1,m*n);
y = reshape(yy,1,m*n);
z = 1+y-.3-3*x;
iz = find(z<0);
niz = length(iz);
z(iz) = zeros(1,niz);
iq = find(y-.3>3*x);
niq = length(iq);
z(iq) = ones(1,niq);
r = (3*x-y-2.2).^2+(x-y-.6).^2;
z = z+ones(1,N)./(1+80*r);
z = reshape(z,m,n);

end