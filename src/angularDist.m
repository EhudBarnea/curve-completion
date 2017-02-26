function angDist = angularDist(or1, or2)
% calculate the angular distance between each value in the vector or1,
% and the scalar or2
angDist = mod(or1 - or2,2*pi);
angDist(angDist > pi) = 2*pi - angDist(angDist > pi);
end
