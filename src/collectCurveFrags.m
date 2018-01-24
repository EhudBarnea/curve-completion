function [] = collectCurveFrags(params)
% Collect curve fragments exhaustively (needs several GB of memory)


% Accumulate all curve fragments into bins relative to a cannonical pose

debug = false;

% minimum distance between two points for which we keep the curve
% fragments.
minDist = -1; % in pixels (Eucledean)

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

% for i = 1:numImgs
parfor i = 1:numImgs
    disp([num2str(i) ' started']);
    
    frags = cell(numBins(1), numBins(2), numOrBins);
    
    % keep the number of samples (edge fragments) for each end point.
    numSamples = zeros(numBins(1), numBins(2), numOrBins);
    
    % go over the different annotators
    for annotatorNum=1:3
        
        % load image curves
        imgName = imgNames{i};
        baseName = imgName(1:end-4);
        % img = imread([imgsFolder imgName]);
        data = load([curvesFolder baseName],'groundTruth');
        curves = data.groundTruth{annotatorNum}; % use one of the annotators
        
        
        % go over all curves.
        % A curve is an Nx2 list of ordered points where the first
        % element is the x coordinate. curves are a list of points, not pixels,
        % so there may be a gap between adjacent points.
        for j=1:length(curves)
            if mod(j,10) == 0
                disp([num2str(i) ',' num2str(annotatorNum) ' - ' num2str(j) '/' num2str(length(curves))]);
            end
            
            % get curve
            c = fixCurve(curves{j}, imgSizes(i,:));
            
            % get curve fragments between each 2 curve points. each point in
            % turn serves as a reference and the curve fragment is organized
            % according to it
            for p1=1:size(c,1)
                for p2=1:size(c,1)
                    if p1==p2
                        continue;
                    end
                    % get curve fragment in canonincal pose
                    [fragPts, s, endPointOr, fragPtsAll, np1, np2, p1o, p2o] = getCanonCurve(c, p1, p2);
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
                    if endPointBin(1) > size(frags,1) || endPointBin(2) > size(frags,2)
                        continue;
                    end
                    
                    % get end point orientation bin
                    endPointOrBin = getOrBin(endPointOr, orBinSize, numOrBins);
                    
                    % place fragment in frags array according to the end point
                    frags{endPointBin(1),endPointBin(2),endPointOrBin} = [frags{endPointBin(1),endPointBin(2),endPointOrBin}; [i, j, p1, p2, endPoint(1), endPoint(2), endPointOr, annotatorNum]];
                    
                    if debug
                        %                     endPoint
                        %                     endPointOrBin
                        
                        close all
                        line(fragPts(:,1),fragPts(:,2));
                        hold on
                        plot(fragPtsAll(p1o:np1,1),fragPtsAll(p1o:np1,2),'r');
                        hold on
                        plot(fragPtsAll(np2:p2o,1),fragPtsAll(np2:p2o,2),'r');
                        axis equal
                        axis([-10 100 -50 50])
                        %                     axis([-200 200 -200 200])
                        saveas(gcf,sprintf('%stmp/%d_%d_%d_%d.png',outFolder,i,j,p1,p2));
                        close all
                    end
                end
            end
        end
    end
    
    % make sure each bin in frags contains only one fragment from each
    % image curve. this is performed after adding all the fragments so that
    % we could choose a fragment randomly from all those that come from the
    % same curve. This can be changed to allow taking several fragments
    % from a single contour if they are far enough from each other.
    for x=1:numBins(1)
        for y=1:numBins(2)
            for z=1:numOrBins
                f = frags{x,y,z};
                if size(f,1) == 0
                    continue;
                end
                diffCurves = unique([f(:,2), f(:,8)],'rows');
                numDiffCurves = size(diffCurves,1);
                filtered = zeros(numDiffCurves,size(f,2));
                for k=1:numDiffCurves
                    % pick random fragment
                    curveFrags = find(f(:,2)==diffCurves(k,1) & f(:,8)==diffCurves(k,2));
                    pickedFrag = f(curveFrags(randi(length(curveFrags))),:);
                    filtered(k,:) = pickedFrag;
                end
                
                frags{x,y,z} = filtered;
                numSamples(x,y,z) = numDiffCurves;
            end
        end
    end
    
    disp([num2str(i) ' saving']);
    parSave([outFolder 'all_frags/frags_from_arr_' num2str(i)], frags, numSamples)
    disp([num2str(i) ' done']);
end


% combine the fragments from all the different images (each in its own
% file) to a single structure
frags = cell(numBins(1), numBins(2), numOrBins);
numSamples = zeros(numBins(1), numBins(2), numOrBins);
for i=1:1:numImgs
    disp(i);
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

disp('saving')
save([outFolder 'all_frags/frags.mat'],'frags','numSamples','-v7.3');
disp('done')

end

