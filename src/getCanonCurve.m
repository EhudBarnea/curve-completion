function [pts] = getCanonCurve(c, p1, p2)

% get curve points pts (in cannonical pose) from curve c relative to curve
% start p1 and until the end p2. p1 can be after p2.

% calculate orientation as the vector connecting a point
% to the next point (or a later point for smoother orientation)

gapSize = 3;
if p1 < size(c,1)
    orVec = c(p1+gapSize,:) - c(p1,:);
else % is last point
    orVec = c(p1,:) - c(p1-1,:);
end
orVec = orVec/norm(orVec);
if p2 < p1
    % flip orientation vector if going from p2 to p1
    orVec = -orVec;
end
% get angle between the orientation vector and the x axis
or = atan2(orVec(2),orVec(1));
or(or<0) = or(or<0) + 2*pi;


% get curve fragment between p1 and p2 relative to p1 as a
% reference. if p2 is before p1 flip the order of curve points
if p1 < p2
    fragPts = c(p1:p2,:); % fragment points
else
    fragPts = c(p1:-1:p2,:); % fragment points (flipped)
end
[fragPts, ~] = transPoints(fragPts, [], c(p1,:), or);
endPoint = fragPts(end,:);

% mirror the curve down.
% if the end point is in the half space above the first
% point (Y>0 relative to the first point), mirror the curve
% to the other side
if endPoint(2) > 0
    fragPts(:,2) = -fragPts(:,2);
%     fragOrs = 2*pi - fragOrs;
    endPoint = fragPts(end,:);
end


pts = fragPts;
end

