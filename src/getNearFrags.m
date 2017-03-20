function [nearFrags] = getNearFrags(endPoint, endPointOr, method, frags, params)
% Get all curve fragments that end 'near' an end point by different methods or
% definitions of 'near'.

% endPoint and endPointOr are assumed to be in cannonical pose

% method:
% method='bin' for fragments with end point inside the same bin.
% method='rad' for fragments with end point up to a certain radius away.
% method='si' for fragments with end point up to a certain radius away and
% in a scale invariant way. For this method supply the scale invariant
% frags dataset in variable 'frags'.

% nearFrags columns:
% 1 = img num
% 2 = curve num
% 3 = start point index
% 4 = end point index

if strcmp(method,'bin')
    nearFrags = getNearFragsBin(endPoint, endPointOr, frags, params);
elseif strcmp(method,'rad')
    nearFrags = getNearFragsRad(endPoint, endPointOr, frags, params);
elseif strcmp(method,'si')
            nearFrags = getNearFragsSI(endPoint, endPointOr, frags, params);
%     nearFrags = getNearFragsSI_withCurva(endPoint, endPointOr, frags, params);
else
    nearFrags = [];
    disp('Error: no such method');
end



end



function nearFrags = getNearFragsBin(endPoint, endPointOr, frags, params)
% Get all curve fragments near an end point. More specifically, those
% for which the endPoint falls in the same bin

if numel(size(frags))~=3
    disp('Error: incorrect frags datastructure');
end

endPointBin = floor((endPoint - [params.relMinX, params.relMinY])/params.binSize) + 1;
endPointBin(endPointBin>params.numBins) = params.numBins(1);
endPointOrBin = getOrBin(endPointOr, params.orBinSize, params.numOrBins);
nearFrags = frags{endPointBin(1), endPointBin(2), endPointOrBin};
end

function nearFrags = getNearFragsRad(endPoint, endPointOr, frags, params)
% Get all curve fragments near an end point. More specifically, those
% that are matchDist away from endPoint and with orientation
% params.matchOr away from endPointOr, where matchDist is calculate
% relative to the distance to the endPoint via the
% params.matchDistFactor parameter.

if numel(size(frags))~=3
    disp('Error: incorrect frags datastructure');
end

% calculate matcing distance (the meaning of "near")
if params.relMatchDist
    inducerDist = norm([0,0] - endPoint);
    matchDist = inducerDist / params.matchDistFactor;
else
    matchDist = params.matchDistFactor;
end

% get the box that bounds the circle of matchDist radius,
% centered at endPoint
boxX = [floor(endPoint(1)-matchDist), ceil(endPoint(1)+matchDist)];
boxY = [floor(endPoint(2)-matchDist), ceil(endPoint(2)+matchDist)];

% move to frags array coordinates
boxX = (boxX - params.relMinX) / params.binSize + 1;
boxY = (boxY - params.relMinY) / params.binSize + 1;

% make sure we don't go out of bounds
boxX(boxX<1) = 1;
boxY(boxY<1) = 1;
boxX(boxX>params.numBins(1)) = params.numBins(1);
boxY(boxY>params.numBins(2)) = params.numBins(2);

% get relevant (somewhat close) frags
relevantFrags = cat(1,frags{boxX(1):boxX(2),boxY(1):boxY(2),:});

% relevantFrags columns:
% 1. Image number
% 2. Curve number
% 3. Start point number
% 4. End point number
% 5. End point X relative to start
% 6. End point Y relative to start
% 7. End point orientation relative to start
% 8. Annotator number

% filter (keep) the relevant (close enough) frags
nearFragsID = normRows(relevantFrags(:,5:6)-repmat(endPoint,size(relevantFrags,1),1)) < matchDist & ...
    angularDist(relevantFrags(:,7), endPointOr) < params.matchOr;
nearFrags = relevantFrags(nearFragsID,:);

% remove fragments from the same curve
[~, idx, ~] = unique(nearFrags(:,[1,2,8]),'rows');
nearFrags = nearFrags(idx,:);
end

function nearFrags = getNearFragsSI(endPoint, endPointOr, frags, params)
% Get all curve fragments near an end point in a scale invariant way.
% More specifically, those
% that are matchDist away from endPoint and with orientation
% params.matchOr away from endPointOr, where matchDist is calculate
% relative to the distance to the endPoint via the
% params.matchDistFactor parameter.

if numel(size(frags))~=2
    disp('Error: incorrect frags datastructure');
end

if params.relMatchDist ~= true
    disp('Error: not implemented for params.relMatchDist=False');
end

% Get angle of vector to endPoint from X axis
angleFromX = 2*pi+atan2(endPoint(2),endPoint(1));
if angleFromX > 2*pi % correct for endPoint=[x<0,0]
    angleFromX = pi;
end

% Get angular limit corresponding to params.matchDistFactor. This is
% the angle of a right triangle with a leg X that is equal to the scale
% and another leg that is equal to X/matchDistFactor.
matchAngle = atan(1/params.matchDistFactor);
matchAngularLimits = [angleFromX-matchAngle, angleFromX+matchAngle];
if matchAngularLimits(1) < pi
    matchAngularLimits(1) = pi;
end
if matchAngularLimits(2) > 2*pi
    matchAngularLimits(2) = 2*pi;
end

% Get bins for the angular limit
matchAngularLimitsBins = floor((matchAngularLimits-pi)/params.siAngBinSize) + 1;
matchAngularLimitsBins(matchAngularLimitsBins > params.siNumAngBins) = params.siNumAngBins;

% Get relevant frags
relevantFrags = cat(1,frags{matchAngularLimitsBins(1):matchAngularLimitsBins(2),:});

% relevantFrags columns:
% 1. Image number
% 2. Curve number
% 3. Start point number
% 4. End point number
% 5. End point X relative to start
% 6. End point Y relative to start
% 7. End point orientation relative to start
% 8. Annotator number
% 9. Angle of vector to endPoint from X axis

% filter (keep) the relevant (close enough) frags
nearFragsID = relevantFrags(:,9) > matchAngularLimits(1) & ...
    relevantFrags(:,9) < matchAngularLimits(2) & ...
    angularDist(relevantFrags(:,7), endPointOr) < params.matchOr;
nearFrags = relevantFrags(nearFragsID,:);

% remove fragments from the same curve
[~, idx, ~] = unique(nearFrags(:,[1,2,8]),'rows');
nearFrags = nearFrags(idx,:);

end


function nearFrags = getNearFragsSI_withCurva(endPoint, endPointOr, frags, params)
% Get all curve fragments near an end point in a scale invariant way.
% More specifically, those
% that are matchDist away from endPoint and with orientation
% params.matchOr away from endPointOr, where matchDist is calculate
% relative to the distance to the endPoint via the
% params.matchDistFactor parameter. Also requires similar curvature at the
% inducer points

% This implementation is incomplete - it assumes that the curvature at the
% end of the to-be-completed curve is 0 (straight line)

nearFrags = getNearFragsSI(endPoint, endPointOr, frags, params);

indCurva = [0,0]; % line
% the arc length to use for the curvature calculation (by circle fitting)
% is determined is the distance between the two inducers divided by this
% amount.
cirlceFitArcLenFactor = 6;
minPtsForFit = 4;

sameCurva = false(size(nearFrags,1),1);
center = zeros(size(nearFrags,1),2);
toUse = size(nearFrags,1);
% toUse = 1000;
for i=1:toUse
    imgNum = nearFrags(i,1);
    cNum = nearFrags(i,2); % curve num
    fragP1 = nearFrags(i,3);
    fragP2 = nearFrags(i,4);
    
    % load curves
    imgName = params.imgNames{imgNum};
    baseName = imgName(1:end-4);
    data = load([params.curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{1}; % use set 1 (each set is a different annotator)
    c = fixCurve(curves{cNum}, params.imgSizes(imgNum,:));
    
    % transform to canonincal pose
    [fragPts,~,~,fragPtsAll] = getCanonCurve(c, fragP1, fragP2);
    
    
    % stretch
    fragPtsAll = [fragPtsAll; fragPts(end,:)];
    fragPtsAll = stretchCurve(fragPtsAll, endPoint);
    fragPtsAll = fragPtsAll(1:end-1,:);
    fragPts = stretchCurve(fragPts, endPoint);
    
    center(i,:) = fragPts(floor(size(fragPts,1)/2),:);
    
    % get curvature points (to calculate the curvature between them)
    cvp1=[0,0];
    cvp2=[0,0];
    if fragP1 < fragP2
        cvp1(2) = fragP1;
        cvp2(1) = fragP2;
    else
        cvp1(2) = size(fragPtsAll,1)-fragP1+1;
        cvp2(1) = size(fragPtsAll,1)-fragP2+1;
    end
    arcLenToP1 = getArcLength(fragPtsAll,cvp1(2),1);
    arcLenFromP2 = getArcLength(fragPtsAll,cvp2(1),size(fragPtsAll,1));
    cirlceFitArcLen = norm(fragPts(end,:))/cirlceFitArcLenFactor;
    cirlceFitArcLen=25;
    cvp1(1) = find(arcLenToP1 < cirlceFitArcLen,1);
    cvp2(2) = find(arcLenFromP2 < cirlceFitArcLen,1,'last');
    
    % make sure we have enough points (3 are enough, but more points = less noise)
    hasCurva = true;
    numPtsSec1 = cvp1(2)-cvp1(1)+1;
    numPtsSec2 = cvp2(2)-cvp2(1)+1;
    if numPtsSec1 < minPtsForFit || numPtsSec2 < minPtsForFit
        hasCurva = false;
    end
    
    % calculate curvature
    if hasCurva
        % fit circle
        [xc1,yc1,R1,~,isLine1] = circfit(fragPtsAll(cvp1(1):cvp1(2),1),fragPtsAll(cvp1(1):cvp1(2),2));
        [xc2,yc2,R2,~,isLine2] = circfit(fragPtsAll(cvp2(1):cvp2(2),1),fragPtsAll(cvp2(1):cvp2(2),2));
        % set curvature (for normal configurations or colinear points)
        if ~isLine1
            curva1 = 1/R1;
        else
            curva1 = 0;
        end
        if ~isLine2
            curva2 = 1/R2;
        else
            curva2 = 0;
        end
        
        % check if curvature is similar
        th = 0.02;
        if abs(curva1-indCurva(1))<th && abs(curva2-indCurva(2))<th
            sameCurva(i) = true;
        end
        
        % ------- vis
        % Visualize curve and curvature
        if ~sameCurva(i) || true
            continue;
        end
        close all
        plot(fragPtsAll(:,1),fragPtsAll(:,2));
        hold on
        scatter(fragPtsAll(cvp1(2),1),fragPtsAll(cvp1(2),2),20)
        hold on
        scatter(fragPtsAll(cvp2(1),1),fragPtsAll(cvp2(1),2),20)
        
        th = linspace(0,2*pi,80)';
        if R1<100
            xe = R1*cos(th)+xc1;
            ye = R1*sin(th)+yc1;
            plot([xe;xe(1)],[ye;ye(1)]);
        end
        hold on
        if R2<100
            xe = R2*cos(th)+xc2;
            ye = R2*sin(th)+yc2;
            plot([xe;xe(1)],[ye;ye(1)]);
        end
        fprintf('%f,%f\n',1/R1,1/R2);
        
        plot(fragPtsAll(cvp1(1):cvp1(2),1),fragPtsAll(cvp1(1):cvp1(2),2),'r');
        hold on
        plot(fragPtsAll(cvp2(1):cvp2(2),1),fragPtsAll(cvp2(1):cvp2(2),2),'r');
        hold on
        axis([-200 200 -200 200]);
        %         axis equal
        
        LogicalStr = {'Circle', 'Line'};
        title(sprintf('%s %d,%s %d',LogicalStr{isLine1+1},numPtsSec1,LogicalStr{isLine2+1},numPtsSec2));
        saveas(gcf,['/Users/ehud/Downloads/curva/tmp/' num2str(i) '.png']);
        close all
        % ------- end vis
    end
end

% ---- vis
close all

subplot(1,2,1)
meanCenter = mean(center,1);
centersSTD = mean(normRows(center-repmat(meanCenter,size(center,1),1)));
plot(center(:,1),center(:,2),'b.');
hold on
plot(center(sameCurva,1),center(sameCurva,2),'g.');
hold on
visInducers([0,0],0,endPoint,endPointOr,false);
hold on
scatter(meanCenter(1),meanCenter(2),30,'r','filled');
hold on;
viscircles(meanCenter,centersSTD(1), 'LineWidth',1);
hold on;
axis([-200 200 -200 200]);
title(sprintf('(%f,%f) %f',meanCenter(1),meanCenter(2),centersSTD));
hold on


subplot(1,2,2)
meanCenter = mean(center(sameCurva,:),1);
centersSTD = mean(normRows(center(sameCurva,:)-repmat(meanCenter,size(center(sameCurva,:),1),1)));
plot(center(sameCurva,1),center(sameCurva,2),'b.');
hold on
visInducers([0,0],0,endPoint,endPointOr,false);
hold on
scatter(meanCenter(1),meanCenter(2),30,'r','filled');
hold on;
viscircles(meanCenter,centersSTD(1), 'LineWidth',1);
hold on;
axis([-200 200 -200 200]);
title(sprintf('(%f,%f) %f',meanCenter(1),meanCenter(2),centersSTD));
% ---- end vis

% filter curves based on curvature
nearFrags = nearFrags(sameCurva,:);

end



function res = normRows(mat)
% calculate the norm of each row in mat
res = sqrt(sum(mat.^2,2));
end

