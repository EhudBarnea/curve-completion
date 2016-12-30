function [or] = getOrFromBin(ob, orBinSize)
% get continuous orientation as the center of orientation bin

% this assume that the first orientation bin is centered at orientation=0

or = mod((ob-1)*orBinSize,2*pi);

end

