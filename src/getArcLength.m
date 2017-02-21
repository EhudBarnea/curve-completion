function [arcLen] = getArcLength(pts, p1, p2)
% For each point between p1 and p2 calculate the arc lengths from point p1.
% p1 and p2 are indices in pts and by default, p1 and p2 are the first and
% last points, respectively.

if nargin < 3
    p1 = 1;
    p2 = size(pts,1);
end

% flip curve if the direction is wrong
flippedCurve = false;
if p1 > p2
    flippedCurve = true;
    pts = pts(end:-1:1,:);
    p1 = size(pts,1) - p1 + 1;
    p2 = size(pts,1) - p2 + 1;
end

arcLen = zeros(size(pts,1),1);
for i=(p1+1):p2
    d = norm(pts(i,:)-pts(i-1,:));
    arcLen(i) = arcLen(i-1) + d;
end

% flip lengths back
if flippedCurve
    arcLen = arcLen(end:-1:1,:);
end

end

