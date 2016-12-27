function orBin = getOrBin(or, orBinSize, numOrientationBins)
% Get the bin number for a given orientation.

% Adjust the bins such that orientation=0 is the center of the first bin.
% This causes orientations 0, pi/2, pi, and 3pi/2 to be exactly the center
% of some bins (where numOrientationBins=8 and in other cases), which makes for nicer looking figures.

startAngle = -orBinSize/2;

orBin = floor((or-startAngle)/orBinSize) + 1;
if orBin > numOrientationBins
    orBin = 1;
end

end


