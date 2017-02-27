function [intAngle] = getLineIntAngle(p2, or2)
% For a pair of inducers calculate the angle of intersection between the
% two lines that "continue" from the inducers' locations. Returns a
% negative number when the lines to not intersect. The inducers are assumed
% to be in cannonical pose (the first inducer is at point (0,0) with
% orientation 0).

if or2 > 0 && or2 <= pi
    % The angle between the two lines and or2, the orientation of the second
    % inducer, are intersecting angles
    intAngle = or2;
    
    % Make sure that the lines intersect at a point (x,0) where x>=0.
    % Calculation:
    % p2+orVec*s = [x,0]
    % p2(2)+orVec(2)*s = 0
    orVec = [cos(or2), sin(or2)];
    
    s = -p2(2)/orVec(2);
    intersectionPoint = p2+orVec*s;
    if intersectionPoint(1) < 0
        intAngle = -1;
    end
else
    % The lines do not intersect
    intAngle = -1;
end


end

