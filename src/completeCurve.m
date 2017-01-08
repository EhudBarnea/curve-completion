function [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, vis)
% complete curve between points p1,p2 with orientations or1,or2

% output:
% c - the completed curve (points)
% isUsable - whether the data can be trusted but a good completion
% out - output struct containing.
% out.numFrags - total number of fragments observed between the two inducers
% out.numDiffImgs - number of images from which the used fragments were taken from
% out.fragCenters - center points of all fragments observed between the two inducers

% maxFragsToUse = 30;
maxFragsToUse = inf;
maxFragsToShow = 10;
numCurveRepPts = 5;

interpolateCurve = false;

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

% get relevant fragments
endPointFrags = getNearFrags(endPoint, endPointOr, frags, params);

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
isUsable = numFragsToUse>=20;

% prepare output struct
out.fragCenters = fragCenters;
out.numDiffImgs = numDiffImgs;
out.numFrags = numFrags;

% calculate mean curve points
%         coeff = pca(allRepPts);
meanPts = mean(allRepPts,1);
meanPts = reshape(meanPts, numCurveRepPts, 2);

% interpolate curve
c = [0,0; meanPts; endPoint];
if interpolateCurve
    c = interpCurve(c, [0;endPointOr]);
end


if vis
    % draw mean curve
    scatter(meanPts(:,1),meanPts(:,2),7,'r','filled');
%     line(meanPts(:,1),meanPts(:,2),'Color','r');
    plot(c(:,1), c(:,2), 'Color', 'r');
    
    hold on
    scatter(0,0,7,'r','filled')
    hold on
    scatter(endPoint(1),endPoint(2),7,'r','filled')
    axis equal
    axis([-200 200 -200 200])
    title(['Num shown curves = ' num2str(min(maxFragsToShow,numFragsToUse)) '   Num used curves = ' num2str(numFragsToUse) '   num Diff Imgs=' num2str(numDiffImgs)])
end

% mirror back curve
if mirrored
    c(:,2) = -c(:,2);
end

% put curve in the coordinate system of p1 and p2
c = transBackPoints(c, p1, or1);

end


function nearFrags = getNearFrags(endPoint, endPointOr, frags, params)
    % Get all curve fragments near an end point. More specifically, those
    % that are params.matchDist away from endPoint and with orientation
    % params.matchOr away from endPointOr.

    % get the box that bounds the circle of params.matchDist radius,
    % centered at endPoint
    boxX = [floor(endPoint(1)-params.matchDist), ceil(endPoint(1)+params.matchDist)];
    boxY = [floor(endPoint(2)-params.matchDist), ceil(endPoint(2)+params.matchDist)];
    
    % move to frags array coordinates
    boxX = (boxX - params.relMinX) / params.binSize + 1;
    boxY = (boxY - params.relMinY) / params.binSize + 1;
    
    % make sure we don't go out of bounds
    boxX(boxX<1) = 1;
    boxY(boxY<1) = 1;
    boxX(boxX>params.numBins(1)) = params.numBins(1);
    boxY(boxY>params.numBins(2)) = params.numBins(2);
    
    % get relevant frags
    relevantFrags = cat(1,frags{boxX(1):boxX(2),boxY(1):boxY(2),:});
    
    % filter the relevant frags
    nearFrags = normRows(relevantFrags(:,5)-endPoint(1)) < params.matchDist & ...
        normRows(relevantFrags(:,6)-endPoint(2)) < params.matchDist & ...
        angularDist(relevantFrags(:,7), endPointOr) < params.matchOr;
    nearFrags = relevantFrags(nearFrags,:);
end


function angDist = angularDist(or1, or2)
    % calculate the angular distance between each value in the vector or1,
    % and the scalar or2
    angDist = mod(or1 - or2,2*pi);
    angDist(angDist > pi) = 2*pi - angDist(angDist > pi);
end

function res = normRows(mat)
    % calculate the norm of each row in mat
    res = sqrt(sum(mat.^2,2));
end

function c = interpCurve(curvePts, startEndOrs)
    % interpolate a smooth curve for the set of points curveRepPts with
    % orientations for the first and last points startEndOrs. This
    % uses a clamped spline, which is a spline for each conditions are
    % given for the end points in the form of their derivatives.
    
    
    % transform curve points into a function by placing the start and end
    % points on the x axis
    refP = [0, 0];
    refOr = getOrTwoPts(curvePts(1,:), curvePts(end,:));
    ors = [startEndOrs(1); zeros(size(curvePts,1)-2,1); startEndOrs(2)];
    [transPts, transOrs] = transPoints(curvePts, ors, refP, refOr);
    
    % fit a clamped spline to the points
    cs = spline(transPts(:,1)',transPts(:,2)'); % without start/end orientations
%     der = 100;
%     cs = spline(transPts(:,1)',[der, transPts(:,2)', -der]);
    xx = linspace(transPts(1,1),transPts(end,1),100); % generate more points
    yy = ppval(cs,xx);
    c = [xx', yy'];
    
    % transform curve back
    [c] = transBackPoints(c, curvePts(1,:), -pi/2);
end