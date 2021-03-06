function   [xc,yc,R,a,isLine] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991,
% Origin:
% http://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit

% Updates by Ehud: we check before fitting whether the data is all along a
% single line (producing a rank deficient matrix). In this case the
% circle's radius is infinity and the curvature is 0.


x=x(:); y=y(:);
aMat = [x y ones(size(x))];
if rank(aMat) > 2
    a=aMat\(-(x.^2+y.^2));
    xc = -.5*a(1);
    yc = -.5*a(2);
    R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
    isLine = false;
else
    R = inf;
    xc = [];
    yc = [];
    a = [];
    isLine = true;
end
