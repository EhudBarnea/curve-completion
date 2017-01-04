function [] = collectCurveFragsRand(params)
% Collect curve fragments randomly
% this part may need to be updated to work with the other parts

% this code doesn't work and needs to be updated

% Accumulate all curve fragments into bins relative to a cannonical pose.
% Sample curve fragments randomly.

% keep the number of samples (edge fragments) for each end point.
% this is indexed by [y,x] for visualization
numSamples = zeros((relMaxY-relMinY)/binSize, (relMaxX-relMinX)/binSize);

numNeededSamples = 200; % per curve bin

tooLongFrags = 0;


while(true)
    % randomly pick image
    i = floor(rand*numImgs + 1);
    
    % load image curves
    imgName = imgNames{i};
    baseName = imgName(1:end-4);
    imgSize = imgSizes(i,:);
    data = load([curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{1}; % use set 1 (each set is probably a different annotator)
    
    
    % randomly pick a curve.
    % A curve is an Nx2 list of ordered points where the first
    % element is the x coordinate. curves are a list of points, not pixels,
    % so there may be a gap between adjacent points.
    j = floor(rand*length(curves)+1);
    c = fixCurve(curves{j});
    
    % randomly pick end points. never pick the last point of a curve, to
    % make the (later) calculation of orientation simpler
    p1 = floor(rand*(size(c,1)-1) + 1);
    p2 = floor(rand*(size(c,1)-1) + 1);
    
    if p1==p2
        continue;
    end
    startPoint = c(p1,:);
    endPoint = c(p2,:);
    % don't use too short curve fragments
    if norm(endPoint - startPoint) < minDist
        continue;
    end
    
    [fragPts] = getCanonCurve(c, p1, p2);
    relEndPoint = fragPts(end,:); % end point relative to the first
    
    % if the end point is in the half space above the first
    % point (Y>0 relative to the first point), mirror the curve
    % to the other side
    if endPoint(2) > 0
        fragPts(:,2) = -fragPts(:,2);
        relEndPoint = fragPts(end,:);
    end
    
    relEndPointBin = floor((relEndPoint - [relMinX, relMinY])/binSize) + 1;
    
    % check if curve fragment is too long (for the array)
    if relEndPoint(1)>relMaxX || relEndPoint(1)<relMinX || relEndPoint(2)>relMaxY || relEndPoint(2)<relMinX
        tooLongFrags = tooLongFrags + 1;
        continue;
    end
    
    % keep number of samples
    flippedY =  size(numSamples,2) - relEndPointBin(2) + 1;
    if numSamples(flippedY, relEndPointBin(1)) < numNeededSamples
        numSamples(flippedY, relEndPointBin(1)) = numSamples(flippedY, relEndPointBin(1)) + 1;
        
        % place fragment in frags array according to the end point
        fragNum = fragNum + 1;
        frags{fragNum,1} = relEndPointBin(1);
        frags{fragNum,2} = relEndPointBin(2);
        frags{fragNum,3} = fragPts;
        
        %         line(fragPts(:,1),fragPts(:,2));
        %         axis equal
        
        display([num2str(relEndPointBin(1)) ' ' num2str(relEndPointBin(1)) ' = ' num2str(numSamples(flippedY, relEndPointBin(1)))]);
    end
end



end

