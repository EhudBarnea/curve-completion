function [pts, s, lastOr, ptsAll] = getCanonCurve(c, p1, p2)

% get curve points pts (in cannonical pose) from curve c relative to curve
% start p1 and until the end p2. p1 can be after p2.

% s - whether the calculation succeeded. It is possible that there aren't
% enough points to calculate the orientation.
% lastOr - orientation of last point.
% ptsAll - all curve points (not only between p1 and p2)

% calculate orientation as the vector connecting a point
% to the point that is gapSize points away
gapSize = 3;

% get points in curve fragment. if p2 is before p1 flip the order of points
fragPts = c;
if p1 > p2
    fragPts = c(end:-1:1,:); % fragment points (flipped)
    p1 = size(fragPts,1)-p1+1;
    p2 = size(fragPts,1)-p2+1;
end

% make sure we have enough points to calculate orientations
numPts = p2 - p1 + 1;
if numPts < gapSize + 1
    s = false;
    pts = [];
    lastOr = 0;
    return;
end


% get orientation p1
firstOr = getOrTwoPts(fragPts(p1,:), fragPts(p1+gapSize,:));

% transform points to be relative to the frame of the first point
[fragPts, ~] = transPoints(fragPts, [], fragPts(p1,:), firstOr);

% mirror the curve down.
% if the end point is in the half space above the first
% point (Y>0 relative to the first point), mirror the curve
% to the other side
lastPoint = fragPts(p2,:);
if lastPoint(2) > 0
    fragPts(:,2) = -fragPts(:,2);
    lastPoint = fragPts(p2,:);
end

% get orientation of last point
lastOr = getOrTwoPts(fragPts(p2,:), fragPts(p2-gapSize,:));

ptsAll = fragPts;
pts = fragPts(p1:p2,:);
s = true;
end

