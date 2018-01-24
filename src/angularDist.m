function angDist = angularDist(or1, or2)
% calculate the angular distance between each value in the scalar/vector or1,
% and the scalar/vector or2
angDist = mod(or1 - or2,2*pi);
angDist(angDist > pi) = 2*pi - angDist(angDist > pi);

% also:
% angDist = acos(cos(or1)*cos(or2)+sin(or1)*sin(or2));
end
