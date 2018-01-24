function [diff, amount] = isDiff(p1, or1, p2, or2)
% Determine whether configuration is difficult

% vector connection p1 to p2
vp1p2 = p2 - p1;
vp1p2 = vp1p2/norm(vp1p2);

% calculate angle1/theta1
v1 = [cos(or1), sin(or1)];
angle1 = acos(v1 * vp1p2');

% calculate angle2/theta2
v2 = [cos(or2), sin(or2)];
angle2 = acos(v2 * vp1p2');

% [angle1, angle2]*180/pi

% determine if difficult
if angle1>pi/2 && angle2<pi/2
    diff = true;
else
    diff = false;
end

% calculate amount of difficulty
d = angularDist([angle1, angle2], [0, pi]);
amount = sum(d);

end

