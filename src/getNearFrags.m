function [nearFrags] = getNearFrags(endPoint, endPointOr, method, frags, params)
% Get all curve fragments that end 'near' an end point by different methods or
% definitions of 'near'.

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
    nearFrags = getNearFragsBin(endPoint, endPointOr, frags, params);
else
    nearFrags = [];
    disp('Error: no such method');
end

end



function nearFrags = getNearFragsBin(endPoint, endPointOr, frags, params)
    % Get all curve fragments near an end point. More specifically, those
    % for which the endPoint falls in the same bin

    
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
    
    % filter (keep) the relevant (close enough) frags
    nearFragsID = normRows(relevantFrags(:,5:6)-repmat(endPoint,size(relevantFrags,1),1)) < matchDist & ...
        angularDist(relevantFrags(:,7), endPointOr) < params.matchOr;
    nearFrags = relevantFrags(nearFragsID,:);
    
    % remove fragments from the same curve
    [~, idx, ~] = unique(nearFrags(:,1:2),'rows');
    nearFrags = nearFrags(idx,:);
end