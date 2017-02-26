function [or] = getOrFit(pts)
% get orientation of the points in pts, pointing to the last point
% pts(end,:), by fitting a line

% fit line to points (unless they form a vertical line)
if any(abs(pts(:,1)-pts(1,1))>0.001)
    p=polyfit(pts(:,1),pts(:,2),1);
    lineFitOr = atan(p(1)); % signed angle from X axis
    if lineFitOr < 0
        lineFitOr = 2*pi + lineFitOr;
    end
else
    lineFitOr = pi/2;
end

% get orientation of the first inducer
orVec = pts(end,:) - pts(1,:);
orVec = orVec/norm(orVec);
% get angle between the orientation vector and the x axis
or = atan2(orVec(2),orVec(1));
if or<0
    or = or + 2*pi;
end

% fix angle of fitted line
if angularDist(lineFitOr,or) > pi/2
    lineFitOr = mod(pi + lineFitOr,2*pi);
end

or = lineFitOr;
end

%%
% 
% pts = [
%     1,1;
%     1,2;
%     1,3;
% 
% ];
% p=polyfit(pts(:,1),pts(:,2),1);
% or = atan(p(1));
% rad2deg(or)