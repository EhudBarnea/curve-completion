function [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, vis)
% complete curve between points p1,p2 with orientations or1,or2

% output:
% c - the completed curve (points)
% isUsable - whether the data can be trusted but a good completion
% out - output struct containing.
% out.numFrags - total number of fragments observed between the two inducers
% out.numDiffImgs - number of images from which the used fragments were taken from
% out.fragCenters - center points of all fragments observed between the two inducers

% maxFragsToUse = 100;
maxFragsToUse = inf;
maxFragsToShow = 10;
numCurveRepPts = 5;

c = [];
out = [];

% vis - visualize completion process
if vis
    figure
end

% get p2 and or2 in relative to p1 and or1
[endPoint, endPointOr] = transPoints(p2, or2, p1, or1);
% mirror p2
mirrored = false;
if endPoint(2) > 0
    mirrored = true;
    endPoint(2) = -endPoint(2);
    endPointOr = 2*pi - endPointOr;
end

endPointOrBin = getOrBin(endPointOr, params.orBinSize, params.numOrBins);

% get relevant fragments
endPointBin = floor((endPoint - [params.relMinX, params.relMinY])/params.binSize) + 1;
endPointBin(endPointBin>params.numBins) = params.numBins(1);
endPointFrags = frags{endPointBin(1), endPointBin(2), endPointOrBin};

numFrags = size(endPointFrags,1);
numFragsToUse = min(numFrags, maxFragsToUse);
if numFrags < 1
    isUsable = false;
    return;
end

% shuffle curves
idx = randperm(numFrags);
endPointFrags = endPointFrags(idx,:);

allRepPts = zeros(numFragsToUse,numCurveRepPts*2); % all curves' representative points
% times 2 because each point is x,y

fragCenters = zeros(numFragsToUse, 2); % center points of all fragments

fragImgs = false(params.numImgs,1); % images with such curves
for i=1:numFragsToUse
    imgNum = endPointFrags(i,1);
    cNum = endPointFrags(i,2); % curve num
    fragP1 = endPointFrags(i,3);
    fragP2 = endPointFrags(i,4);
    
    % load curves
    imgName = params.imgNames{imgNum};
    baseName = imgName(1:end-4);
    data = load([params.curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{1}; % use set 1 (each set is a different annotator)
    c = fixCurve(curves{cNum}, params.imgSizes(imgNum,:));
    
    % transform to canonincal pose
    [fragPts] = getCanonCurve(c, fragP1, fragP2);
    
    fragImgs(imgNum) = true;
    
    % get representative points
    repPts = getCurveEquiPoints(fragPts, numCurveRepPts);
    allRepPts(i,:) = reshape(repPts,1,numCurveRepPts*2);
    
    % get frag center (for further statistics)
    fragCenters(i,:) = getCurveEquiPoints(fragPts, 1);
    
    % display
    if vis && i<=maxFragsToShow
        line(fragPts(:,1),fragPts(:,2));
        hold on
    end
end

% number of different images with such seen curves
numDiffImgs = sum(fragImgs);
% isUsable = numFragsToUse>=20 && numDiffImgs>=5;
isUsable = numFragsToUse>=30;

% prepare output struct
out.fragCenters = fragCenters;
out.numDiffImgs = numDiffImgs;
out.numFrags = numFrags;
out.endPointBin = endPointBin;
out.endPointOrBin = endPointOrBin;

% calculate mean curve
%         coeff = pca(allRepPts);
meanPts = mean(allRepPts,1);
meanPts = reshape(meanPts, numCurveRepPts, 2);

if vis
    % draw mean curve
    scatter(meanPts(:,1),meanPts(:,2),7,'r','filled');
    line(meanPts(:,1),meanPts(:,2),'Color','r');
    
    hold on
    scatter(0,0,7,'r','filled')
    hold on
    scatter(endPoint(1),endPoint(2),7,'r','filled')
    axis equal
    axis([-200 200 -200 200])
    title(['Num shown curves = ' num2str(min(maxFragsToShow,numFragsToUse)) '   num Diff Imgs=' num2str(numDiffImgs)])
end

% mirror back curve
if mirrored
    meanPts(:,2) = -meanPts(:,2);
end

% put curve in the coordinate system of p1 and p2
[meanPts] = transBackPoints(meanPts, p1, or1);
c = [p1; meanPts; p2];

% numFrags
end

