function [] = collectCurveFrags(params)
% Collect curve fragments exhaustively (needs several GB of memory)


% Accumulate all curve fragments into bins relative to a cannonical pose

debug = false;

% minimum distance between two points for which we keep the curve
% fragments.
minDist = 3; % in pixels (Eucledean)



relMinX = params.relMinX;
relMaxX = params.relMaxX;
relMinY = params.relMinY;
relMaxY = params.relMaxY;
numOrBins = params.numOrBins;
binSize = params.binSize;
orBinSize = params.orBinSize;
numBins = params.numBins;
numImgs = params.numImgs;
imgNames = params.imgNames;
imgSizes = params.imgSizes;
curvesFolder = params.curvesFolder;
outFolder = params.outFolder;



parpool(32); % sge59


% for i = 1:1:numImgs
parfor i = 1:numImgs
    display([num2str(i) ' started']);
    
    frags = cell(numBins(1), numBins(2), numOrBins);
    
    % keep the number of samples (edge fragments) for each end point.
    % this is indexed by [y,x] for visualization
    numSamples = zeros(numBins(1), numBins(2), numOrBins);
    
    % load image curves
    imgName = imgNames{i};
    baseName = imgName(1:end-4);
    % img = imread([imgsFolder imgName]);
    data = load([curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{1}; % use set 1 (each set is probably a different annotator)
    
    
    % go over all curves.
    % A curve is an Nx2 list of ordered points where the first
    % element is the x coordinate. curves are a list of points, not pixels,
    % so there may be a gap between adjacent points.
    for j=1:length(curves)
        if mod(j,10) == 0
            display([num2str(i) ' - ' num2str(j) '/' num2str(length(curves))]);
        end
        
        % get curve
        c = fixCurve(curves{j}, imgSizes(i,:));
        
        % get curve fragments between each 2 curve points. each point in
        % turn serves as a reference and the curve fragment is organized
        % according to it
        for p1=1:size(c,1)
            for p2=1:size(c,1)
                % get curve fragment in canonincal pose
                [fragPts, s, endPointOr] = getCanonCurve(c, p1, p2);
                % check if was able to get cannonical curve
                if ~s 
                    continue;
                end
                endPoint = fragPts(end,:);
                
                % check if curve fragment is too short
                if norm(endPoint) < minDist
                    continue;
                end
                
                % check if end point is off limits
                if endPoint(1)>=relMaxX || endPoint(1)<=relMinX || endPoint(2)>=relMaxY || endPoint(2)<=relMinX
                    continue;
                end
                
                % get end point bin
                endPointBin = floor((endPoint - [relMinX, relMinY])/binSize) + 1;
                if endPointBin(1) > size(numSamples,2) || endPointBin(2) > size(numSamples,1)
                    continue;
                end
                
                % get end point orientation bin
                endPointOrBin = getOrBin(endPointOr, orBinSize, numOrBins);
                
                % keep number of samples (in image coordinates (row/col) for
                % visualization)
                flippedY =  size(numSamples,2) - endPointBin(2) + 1;
                numSamples(flippedY, endPointBin(1), endPointOrBin) = numSamples(flippedY, endPointBin(1), endPointOrBin) + 1;
                
                % place fragment in frags array according to the end point
                frags{endPointBin(1),endPointBin(2),endPointOrBin} = [frags{endPointBin(1),endPointBin(2),endPointOrBin}; [i, j, p1, p2]];
                
                if debug
                    endPoint
                    endPointOrBin
                    
                    close all
                    line(fragPts(:,1),fragPts(:,2));
                    axis equal
                    axis([0 100 -50 50])
%                     axis([-200 200 -200 200])
%                     export_fig([outFolder 'tmp/' num2str([i,j,p1,p2])])
%                     close all
                end
            end
        end

    end
    
    display([num2str(i) ' saving']);
    parSave([outFolder 'all_frags/frags_from_arr_' num2str(i)], frags, numSamples)
    display([num2str(i) ' done']);
end



frags = cell(numBins(1), numBins(2), numOrBins);
numSamples = zeros(numBins(1), numBins(2), numOrBins);
for i=1:1:numImgs
    display(i);
    data=load([outFolder 'all_frags/frags_from_arr_' num2str(i)],'frags','numSamples');
    imgFrags = data.frags;
    imgNumSamples = data.numSamples;
    for x=1:numBins(1)
        for y=1:numBins(2)
            for z=1:numOrBins
                frags{x,y,z} = [frags{x,y,z}; imgFrags{x,y,z}];
                numSamples(x,y,z) = numSamples(x,y,z) + imgNumSamples(x,y,z);
            end
        end
    end
end

save([outFolder 'all_frags/frags'],'frags','numSamples','-v7.3');
display('done')

end

