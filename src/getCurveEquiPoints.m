function [pts] = getCurveEquiPoints(curvePts, numPts)
% get a set of numPts points from curvePts such that the arc length
% between them are the same and they are all between the start and end
% points.


% calculate arc length from the beginning to every point
arcLen = getArcLength(curvePts);
curveArcLen = arcLen(end);

% return the point closest to each part of the curve's arc length
pts = zeros(numPts,2);
for i=1:numPts
    arcLenFromStart = i*curveArcLen/(numPts+1);
    [~,pIdx] = min(abs(arcLen - arcLenFromStart));
    p = curvePts(pIdx,:);
    pts(i,:) = p;
end

end

