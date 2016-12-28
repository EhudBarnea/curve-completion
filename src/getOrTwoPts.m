function [or] = getOrTwoPts(p1, p2)
% get orientation from p1 to p2

% get orientation of the first inducer
orVec = p2 - p1;
orVec = orVec/norm(orVec);
% get angle between the orientation vector and the x axis
or = atan2(orVec(2),orVec(1));
if or<0
    or = or + 2*pi;
end

end

