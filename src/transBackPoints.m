function [transPts] = transBackPoints(pts, refP, refOr)
% transform back the set of points pts according to the reference point refP and its orientation refOr

% points are assumed to be relative to point (0,0) already


% amount to rotate (can be negative), so that the reference
% point that looks right will look to direction refOr. (this explanation needs to be made more general for the function interpCurve in completeCurve.m)
rotation_angle = refOr;

% create rotation matrix
R = [cos(rotation_angle), -sin(rotation_angle); sin(rotation_angle) cos(rotation_angle)];
% rotate points
transPts = (R*pts')';

% place points after refP
transPts = transPts + repmat(refP,size(transPts,1),1);

end

