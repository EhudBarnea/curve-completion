
function [c] = fixCurve(c, imgSize)
% fix CFGD dataset curves

% keep only xy coordinates and translate to matlab
% 1-base format.
c = c(:,1:2);
c = c + 1;

% flip y axis (so y axis points up)
c(:,2) = imgSize(2) - c(:,2) + 1;

% remove adjacent equal points (this sometimes happens)
samePts = [all(c(1:end-1,:)==c(2:end,:),2); 0];
c = c(~samePts,:);

end

