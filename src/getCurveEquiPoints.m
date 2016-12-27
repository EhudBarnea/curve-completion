function [pts] = getCurveEquiPoints(curvePts, numPts)
% get a set of numPts points from curvePts such that the geodesic distances
% between them are the same and they are all between the start and end
% points.

pts = zeros(numPts,2);

% calculate geodesic distance from the beginning to every point
geoDist = zeros(size(curvePts,1),1);
for i=2:size(curvePts,1)
    d = norm(curvePts(i,:)-curvePts(i-1,:));
    geoDist(i) = geoDist(i-1) + d;
end
arcLength = geoDist(end);

% return the point closest to each part of the arcLength
for i=1:numPts
    geoDistFromStart = i*arcLength/(numPts+1);
    [~,pIdx] = min(abs(geoDist - geoDistFromStart));
    p = curvePts(pIdx,:);
    pts(i,:) = p;
end

end

