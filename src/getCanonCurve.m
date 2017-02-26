function [pts, s, lastOr, ptsAll, np1, np2, p1o, p2o] = getCanonCurve(c, p1, p2)

% get curve points pts (in cannonical pose) from curve c relative to curve
% start p1 and until the end p2. p1 can be after p2.

% Output:
% s - whether the calculation succeeded. It is possible that there aren't
% enough points to calculate the orientation.
% lastOr - orientation of last point.
% ptsAll - all curve points (not only between p1 and p2)
% np1,np2 - updated p1 and p2 (in cases where p1>p2 and the curve was flipped)
% p1o,p2o - indices of the points used to calculate the orientations at p1 and p2


% calculate orientation by fitting a line using points at most gapSize away
gapSize = 3;


% get points in curve fragment. if p2 is before p1 flip the order of points
fragPts = c;
if p1 > p2
    fragPts = c(end:-1:1,:); % fragment points (flipped)
    p1 = size(fragPts,1)-p1+1;
    p2 = size(fragPts,1)-p2+1;
end

arcLenToP1 = getArcLength(fragPts,p1,1);
arcLenFromP2 = getArcLength(fragPts,p2,size(fragPts,1));
p1o = find(arcLenToP1 < gapSize,1);
p2o = find(arcLenFromP2 < gapSize,1,'last');

% make sure we have enough points to calculate orientations
if p1-p1o<3 || p2o-p2<3
    s = false;
    pts = [];
    ptsAll = [];
    lastOr = 0;
    np1 = p1;
    np2 = p2;
    return;
end

% get orientation p1
firstOr = getOrFit(fragPts(p1o:p1,:));

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
lastOr = getOrFit(fragPts(p2o:-1:p2,:));

ptsAll = fragPts;
pts = fragPts(p1:p2,:);
s = true;
np1 = p1;
np2 = p2;
end

