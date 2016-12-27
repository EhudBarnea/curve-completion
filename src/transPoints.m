function [transPts, transOrs] = transPoints(pts, ors, refP, refOr)
% transform the set of points pts and their orientations ors according to the reference point refP and its orientation refOr

% translate points such that refP is the origin
pts = pts - repmat(refP,size(pts,1),1);

% amount to rotate (can be negative), so that the reference
% point looks up
% rotation_angle = pi/2 - refOr;

% amount to rotate (can be negative), so that the reference
% point looks right
rotation_angle = -refOr;

% create rotation matrix
R = [cos(rotation_angle), -sin(rotation_angle); sin(rotation_angle) cos(rotation_angle)];
% rotate points
transPts = (R*pts')';
% update orientations
transOrs = mod(ors - refOr,2*pi);


end

