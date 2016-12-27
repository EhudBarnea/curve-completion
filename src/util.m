  

%% params


addpath('export_fig/');

datasetFolder = '../data/curve fragments dataset/CFGD_release/';
curvesFolder = [datasetFolder 'GT_mat_CFGD_format/'];
imgsFolder = [datasetFolder 'img/'];

outFolder = '../data/';

% minimum distance between two points for which we keep the curve
% fragments.
minDist = 3; % in pixels (Eucledean)

binSize = 1; % in pixels
% size relative to reference point (adjusted to be a multiple of binSize)
relMaxX = 200;
relMaxY = 200;
relMaxX = floor(relMaxX/binSize)*binSize;
relMaxY = floor(relMaxY/binSize)*binSize;
relMinX = -relMaxX;
relMinY = -relMaxY;
numBins = [(relMaxX-relMinX)/binSize, (relMaxY-relMinY)/binSize];

% get all image names and sizes
files = dir([imgsFolder '*jpg']);
numImgs = length(files);
imgNames = cell(size(files,1),1);
imgSizes = zeros(size(files,1),2);
for i=1:size(files,1)
    imgNames{i} = files(i).name;
    img = imread([imgsFolder files(i).name]);
    imgSizes(i,:) = [size(img,1), size(img,2)];
end

% discretization of end point orientation
numOrBins = 8;
orBinSize = 2*pi/numOrBins;


%% Collect curve fragments exhaustively (needs several GB of memory)


% Accumulate all curve fragments into bins relative to a cannonical pose



% minimum distance between two points for which we keep the curve
% fragments.
minDist = 3; % in pixels (Eucledean)


parpool(32); % sge59


% for i = 1:numImgs
parfor i = 1:numImgs
    
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
    
    % duplicate curves - one from start to end and the other from the
    % end to the start
    numCurves = length(curves);
    for j=1:numCurves
        curves{numCurves+j} = curves{j}(end:-1:1,:);
    end
    numCurves = length(curves);
    
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
        
        % curve number before duplicating curves
        origCurveNum = j;
        if j > numCurves/2
            origCurveNum = j - numCurves/2;
        end
        
        % calculate orientation at each point as the vector connecting it
        % to the point gapSize away from it. the last points has the same orientation as the
        % last point for which we could calculate an orientation.
        gapSize = 3; % in number of points
        if size(c,1) <= gapSize
            continue;
        end
        orVecs = c(1+gapSize:end,:) - c(1:end-gapSize,:);
        orVecs = [orVecs; repmat(orVecs(end,:),gapSize,1)];
        % normalize rows
        rowNorms = sqrt(sum(orVecs.^2,2));
        orVecs = orVecs./[rowNorms, rowNorms];
        % get angle between the orientation vectors and the x axis
        ors = atan2(orVecs(:,2),orVecs(:,1));
        ors(ors<0) = ors(ors<0) + 2*pi;
        
        % get curve fragments between each 2 curve points. each point in
        % turn serves as a reference and the curve fragment is organized
        % according to it. go from beginning to end, then do the opposite
        for p1=1:size(c,1)-1
            for p2=p1+1:size(c,1)
                fragPts = c(p1+1:p2,:); % fragment points
                fragOrs = ors(p1+1:p2);
                [fragPts, fragOrs] = transPoints(fragPts, fragOrs, c(p1,:), ors(p1));
                endPoint = fragPts(end,:);
                
                % check if curve fragment is too short
                if norm(endPoint) < minDist
                    continue;
                end
                
                % check if curve fragment is too long
                if endPoint(1)>=relMaxX || endPoint(1)<=relMinX || endPoint(2)>=relMaxY || endPoint(2)<=relMinX
                    continue;
                end
                
                % mirror the curve down.
                % if the end point is in the half space above the first
                % point (Y>0 relative to the first point), mirror the curve
                % to the other side
                if endPoint(2) > 0
                    fragPts(:,2) = -fragPts(:,2);
                    fragOrs = 2*pi - fragOrs;
                    endPoint = fragPts(end,:);
                end
                
                % get end point bin
                endPointBin = floor((endPoint - [relMinX, relMinY])/binSize) + 1;
                if endPointBin(1) > size(numSamples,2) || endPointBin(2) > size(numSamples,1)
                    continue;
                end
                
                % get end point orientation
                endPointOr = fragOrs(end);
                orBin = getOrBin(endPointOr, orBinSize, numOrBins);
                
                % keep number of samples (in image coordinates (row/col) for
                % visualization)
                flippedY =  size(numSamples,2) - endPointBin(2) + 1;
                numSamples(flippedY, endPointBin(1), orBin) = numSamples(flippedY, endPointBin(1), orBin) + 1;
                
                % place fragment in frags array according to the end point
                frags{endPointBin(1),endPointBin(2),orBin} = [frags{endPointBin(1),endPointBin(2),orBin}; [i, origCurveNum, p1, p2]];
                
%                 orBin
%                 close all
%                 line(fragPts(:,1),fragPts(:,2));
%                 axis equal
%                 axis([0 100 -50 50])

            end
        end

    end
    
    display([num2str(i) ' saving']);
    parSave([outFolder 'all_frags/frags_from_arr_' num2str(i)], frags, numSamples)
    display([num2str(i) ' done']);
end


frags = cell(numBins(1), numBins(2), numOrBins);
numSamples = zeros(numBins(1), numBins(2), numOrBins);
for i=1:numImgs
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


%% Collect curve fragments randomly
% this part may need to be updated to work with the other parts

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

%%
imagesc(numSamples)
axis equal
figure
imagesc(numSamples>0)
axis equal
tooLongFrags

%% For each two inducers show all curves and mean curve

numCurveRepPts = 5;
vis = true;
minDist = 3;

if ~exist('frags')
    load([outFolder 'all_frags/frags']);
end

for x=10%relMinX:10:relMaxX
    for y=-70%relMinY:10:0 % end at 0, because we mirror curves in y>0
        [x, y]
        
        % ignore too short curves (there's too many of them)
        if x<=minDist && x>=-minDist && y<=minDist && y>=-minDist
            continue;
        end
        
        for ob=3%1:numOrBins
            endPoint = [x,y];
            endPointBin = floor((endPoint - [relMinX, relMinY])/binSize) + 1;
            endPointBin(endPointBin>numBins) = numBins(1);
            endPointFrags = frags{endPointBin(1), endPointBin(2), ob};
            
            numFrags = size(endPointFrags,1);
            if numFrags < 1
                continue;
            end
            
            % shuffle curves
            idx = randperm(numFrags);
            endPointFrags = endPointFrags(idx,:);
            
            allRepPts = zeros(numFrags,numCurveRepPts*2); % all curves' representative points
            % times 2 because each point is x,y
            
            maxFragsToShow = 40;
            fragImgs = false(numImgs,1); % images with such curves
            for i=1:numFrags
                imgNum = endPointFrags(i,1);
                cNum = endPointFrags(i,2); % curve num
                p1 = endPointFrags(i,3);
                p2 = endPointFrags(i,4);
                
                % load curves
                imgName = imgNames{imgNum};
                baseName = imgName(1:end-4);
                data = load([curvesFolder baseName],'groundTruth');
                curves = data.groundTruth{1}; % use set 1 (each set is a different annotator)
                c = fixCurve(curves{cNum}, imgSizes(imgNum,:));
                
                % transform to canonincal pose
                [fragPts] = getCanonCurve(c, p1, p2);
                
                fragImgs(imgNum) = true;
                
                % get representative points
                repPts = getCurveEquiPoints(fragPts, numCurveRepPts);
                allRepPts(i,:) = reshape(repPts,1,numCurveRepPts*2);
                
                % display
                if vis && i<maxFragsToShow
                    line(fragPts(:,1),fragPts(:,2));
                    hold on
                end
            end
            
            % number of different images with such seen curves
            numDiffImgs = sum(fragImgs);
            
            % calculate mean curve
            %         coeff = pca(allRepPts);
            meanPts = mean(allRepPts,1);
            meanPts = reshape(meanPts, numCurveRepPts, 2);
            meanPts = [[0,0]; meanPts; endPoint];
            
            meanCurve = meanPts;
%             save(['curves_between_inducers/' num2str(endPoint(1)) '_' num2str(endPoint(2)) '_' num2str(orBin) '.mat'],'meanCurve');
            
            if vis && numDiffImgs>1
                % draw mean curve
                scatter(meanPts(:,1),meanPts(:,2),7,'r','filled');
                line(meanPts(:,1),meanPts(:,2),'Color','r');
                
                hold on
                scatter(0,0,7,'r','filled')
                hold on
                scatter(endPoint(1),endPoint(2),7,'r','filled')
                axis equal
                axis([-200 200 -200 200])
                title(['Num shown curves = ' num2str(min(40,numFrags)) '   num Diff Imgs=' num2str(numDiffImgs)])
                
                %export_fig([outFolder '/curves_between_inducers/c_' num2str(endPoint(1)) '_' num2str(endPoint(2)) '_' num2str(ob) '.png']);
            end
            %close all;
        end
    end
end

