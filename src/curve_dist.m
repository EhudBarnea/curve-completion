function [d, dRel] = curve_dist(c1, c2)
% calculate absolute and relative distance between curves.
% dRel is the distance relative to the arc-lengths of c1.

% distance between the start and end of c1
dC1 = sqrt(sum((c1(end,:) - c1(1,:)) .^ 2)); % Eucledean distance between end points
% dC1 = getArcLength(c1);
% dC1 = dC1(end);

% flip curve c2 if they are in a different order - if the first point of c1
% is closest to the first or last in c2
d1 = sqrt(sum((c1(1,:) - c2(1,:)) .^ 2));
d2 = sqrt(sum((c1(1,:) - c2(end,:)) .^ 2));
if d2 < d1
    c2 = c2(size(c2,1):-1:1,:);
end

% calculate distance between curves
% d = maxPointPairDist(c1, c2); % maximal distance
[d, ~] = DiscreteFrechetDist(c1, c2); % Frechet distance
dRel = d / dC1;
end

function [d] = maxPointPairDist(c1, c2)
% calculate distance d between curves as the maximal distance between closest
% points in c1 and c2 (calculated both from c1 to c2 and from c2 to c1).

% find closest points in gt and in c (naive/slow implementation)
c1Closest = zeros(size(c1,1),1); % index of closest c2 point
c1ClosestD = inf(size(c1,1),1); % distance to closest c2 point
c2Closest = zeros(size(c2,1),1);
c2ClosestD = inf(size(c2,1),1);
for j=1:size(c1,1)
    for k=1:size(c2,1)
        d = sqrt(sum((c1(j,:) - c2(k,:)) .^ 2));
        if d < c1ClosestD(j)
            c1ClosestD(j) = d;
            c1Closest(j) = k;
        end
        if d < c2ClosestD(k)
            c2ClosestD(k) = d;
            c2Closest(k) = j;
        end
    end
end

% get maximal distance between closest points
d = max([c1ClosestD; c2ClosestD]);

end
