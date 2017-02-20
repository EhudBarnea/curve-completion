function pts = stretchCurve(pts, endPoint)
% Stretch (and rotate) curve points pts such that the last point in pts will match
% endPoint

% Rotate pts such that pts(end,:) is on the X axis
theta = -atan2(pts(end,2),pts(end,1));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pts = (R*(pts'))';

% Stretch pts in both directions according to the distance to endPoint
stretchFactor = norm(endPoint)/norm(pts(end,:));
pts = pts.*repmat(stretchFactor,size(pts,1),2);

% Put curve back in place according to the original location of endPoint
theta = atan2(endPoint(2),endPoint(1));
rad2deg(theta);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pts = (R*(pts'))';

end
