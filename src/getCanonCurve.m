function [pts, s, lastOr] = getCanonCurve(c, p1, p2)

% get curve points pts (in cannonical pose) from curve c relative to curve
% start p1 and until the end p2. p1 can be after p2.

% s - whether the calculation succeeded. It is possible that there aren't
% enough points to calculate the orientation.

% lastOr - orientation of last point

% calculate orientation as the vector connecting a point
% to the point that is gapSize points away
gapSize = 3;

% get points in curve fragment. if p2 is before p1 flip the order of points
if p1 < p2
    fragPts = c(p1:p2,:); % fragment points
else
    fragPts = c(p1:-1:p2,:); % fragment points (flipped)
end

% make sure we have enough points to calculate orientations
numPts = size(fragPts,1);
if numPts < gapSize + 1
    s = false;
    pts = [];
    lastOr = 0;
    return;
end


% get orientation of the first point
orVec = fragPts(1+gapSize,:) - fragPts(1,:);
orVec = orVec/norm(orVec);
% get angle between the orientation vector and the x axis
or = atan2(orVec(2),orVec(1));
if or<0
    or = or + 2*pi;
end

% transform points to be relative to the frame of the first point
[fragPts, ~] = transPoints(fragPts, [], fragPts(1,:), or);

% mirror the curve down.
% if the end point is in the half space above the first
% point (Y>0 relative to the first point), mirror the curve
% to the other side
lastPoint = fragPts(end,:);
if lastPoint(2) > 0
    fragPts(:,2) = -fragPts(:,2);
    lastPoint = fragPts(end,:);
end

% get orientation of last point
orVec = fragPts(end-gapSize,:) - fragPts(end,:);
orVec = orVec/norm(orVec);
% get angle between the orientation vector and the x axis
or = atan2(orVec(2),orVec(1));
if or<0
    or = or + 2*pi;
end
lastOr = or;

pts = fragPts;
s = true;
end

