function [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, vis)
% complete curve between points p1,p2 with orientations or1,or2

% input:
% frags - frags dataset. use fragsSI if params.matchSI=true
% output:
% c - the completed curve (points)
% isUsable - whether the data can be trusted but a good completion
% out - output struct containing.
% out.numFrags - total number of fragments observed between the two inducers
% out.numDiffImgs - number of images from which the used fragments were taken from
% out.fragRepPts - rep points of all fragments observed between the two inducers
% out.fragCenters - center points of all fragments observed between the two inducers

% maxFragsToShow = 10;
maxFragsToShow = 5;


c = [];
out = [];

calcCenters = false;
if nargout > 3
    calcCenters = true;
end
    

% vis - visualize completion process
if vis
%     figure
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
if ~params.matchSI
    endPointFrags = getNearFrags(endPoint, endPointOr, 'rad', frags, params);
else
    endPointFrags = getNearFrags(endPoint, endPointOr, 'si', frags, params);
end

% remove frags from small scales
if params.matchSI
    fragScales = normRows(endPointFrags(:,5:6));
    endPointFrags = endPointFrags(fragScales > params.siMinScale,:);
end

numFrags = size(endPointFrags,1);
numFragsToUse = min(numFrags, params.maxFragsToUse);
if numFrags < 1
    isUsable = false;
    return;
end

% shuffle curves
idx = randperm(numFrags);
endPointFrags = endPointFrags(idx,:);

allRepPts = zeros(numFragsToUse,params.numCurveRepPts*2); % all curves' representative points
% times 2 because each point is x,y

fragCenters = zeros(numFragsToUse, 2); % center points of all fragments

fragImgs = false(params.numImgs,1); % images with such curves
for i=1:numFragsToUse
    imgNum = endPointFrags(i,1);
    cNum = endPointFrags(i,2); % curve num
    fragP1 = endPointFrags(i,3);
    fragP2 = endPointFrags(i,4);
    annotatorNum = endPointFrags(i,8);
    
    % load curves
    imgName = params.imgNames{imgNum};
    baseName = imgName(1:end-4);
    data = load([params.curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{annotatorNum}; % use one of the annotators
    c = fixCurve(curves{cNum}, params.imgSizes(imgNum,:));
    
    % transform to canonincal pose
    [fragPts] = getCanonCurve(c, fragP1, fragP2);
    
    % stretch such that the end point matches endPoint
    if params.relMatchDist || params.matchSI
        fragPts = stretchCurve(fragPts, endPoint);
    end
    
    fragImgs(imgNum) = true;
    
    % get representative points
    repPts = getCurveEquiPoints(fragPts, params.numCurveRepPts);
    allRepPts(i,:) = reshape(repPts,1,params.numCurveRepPts*2);
    
    % get frag center
    fragCenters(i,:) = getCurveEquiPoints(fragPts, 1);
    
    % display
    if vis && i<=maxFragsToShow
%         line(fragPts(:,1),fragPts(:,2),'Color','k');
%         hold on
        % ---- show each fragment in its own image
        %         axis equal
        %         axis([-200 200 -200 200])
        %         pause(0.5)
        %         close all
        % ----
    end
end

% number of different images with such seen curves
numDiffImgs = sum(fragImgs);
isUsable = numFragsToUse>=20;

% prepare output struct
out.fragRepPts = allRepPts;
out.fragCenters = fragCenters;
out.numDiffImgs = numDiffImgs;
out.numFrags = numFrags;

% calculate completion curve
c = genCurveMean(allRepPts, params);
c = [0,0; c; endPoint];

% split the curve into two parts if not enough samples
if params.extensible && numFrags<params.numMinFrags
%     fprintf('calculating extensible curve, %d\n',numFrags);
    % Get mean center points and orientations
    splitPointIdx = floor(size(c,1)/2)+1;
    splitPoint = c(splitPointIdx,:);
    d=10;
    splitPointBefore = c(splitPointIdx-d,:);
    splitPointAfter = c(splitPointIdx+d,:);
    splitOr1 = getOrTwoPts(splitPointAfter, splitPointBefore);
    splitOr2 = getOrTwoPts(splitPointBefore, splitPointAfter);
    
    paramsExt = params;
    paramsExt.extensible = false;
    paramsExt.matchOr = paramsExt.extMatchOr;
    [c1,~,out1] = completeCurve([0,0], 0, splitPoint, splitOr1, frags, paramsExt, false);
    [c2,~,out2] = completeCurve(splitPoint, splitOr2, endPoint, endPointOr, frags, paramsExt, false);
    c = [c1; c2(2:end,:)];
    numFrags = min([out1.numFrags, out2.numFrags]);
%     scatter(splitPoint(1),splitPoint(2),150,'r')
%     hold on
end




if vis
    % draw mean curve
%     scatter(c(:,1),c(:,2),7,'b','filled'); % draw completion points
%     curveColor = [0,0,0.5];
    maxNum = 300;
%     minNum = 30;
    minNum = 0;
    cmap = colormap(winter(maxNum-minNum+1));
    cmap=flipud(cmap);
    colormap(cmap);
    tmpFrags = numFrags;
    if tmpFrags>maxNum
        tmpFrags=maxNum;
    end
    if tmpFrags<minNum
        tmpFrags=minNum;
    end
%     tmpFrags
%     caxis manual
    caxis([minNum maxNum]);
    curveColor = cmap(tmpFrags-minNum+1,:);
    plot(c(:,1), c(:,2), 'Color', curveColor, 'LineWidth' , 1); % draw completed curve
    
    % -- colorbar
    hcb=colorbar;
    tcks = hcb.Ticks;
    tcks=num2cell(tcks);
    tcks{end}=['>' num2str(tcks{end})];
    set(hcb,'TickLabels',tcks)
    % --
    
    hold on
    visInducers([0, 0], 0, endPoint, endPointOr, false);
%     scatter(0,0,7,'r','filled')
%     hold on
%     scatter(endPoint(1),endPoint(2),7,'b','filled')
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

function curvePts = genCurveMean(allRepPts, params)
% Generate curve using a set of curves as their mean curve points, i.e.,
% each point in the curve is a mean of all the other corresponding points.

% allRepPts - the set of representative points of the curves to use. This
% is an array of NxM where N is the number of curves and M is
% params.numCurveRepPts*2 (for x and y).

meanPts = mean(allRepPts,1);
meanPts = reshape(meanPts, params.numCurveRepPts, 2);

curvePts = meanPts;
end



function res = normRows(mat)
% calculate the norm of each row in mat
res = sqrt(sum(mat.^2,2));
end

