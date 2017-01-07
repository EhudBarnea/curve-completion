
% Complete curves in an image (after marking an occluder)

% this code is obsolete and needs to be updated to include orientations and
% other small changes.


imgNum = 13; % car
% imgNum = 4; % kid
% imgNum = 24; % duck
% imgNum = 1; % plane

showEdgeImg = true;


numCurveRepPts = 5;

if ~exist('frags')
    load([params.outFolder 'all_frags/frags']);
end


close all

imgName = imgNames{imgNum};
baseName = imgName(1:end-4);
img = imread([params.imgsFolder imgName]);
imgRGB = img;
% img = rgb2gray(img);
img = img(:,:,1);

% get edges and their orientations
threshold = 0.2;
[edgesBefore,~,edgeOrDX, edgeOrDY] = edgeMod(img,'Canny',threshold); % modified to return edge orientations


% get occluder
imshow(edgesBefore);
set(gcf,'position',get(0,'screensize'))
[x,y] = ginput(2);
close all
box = int32([min(y), min(x), max(y), max(x)]);

% add occluder
% box = [y1,x1,y2,x2]
img(box(1):box(3),box(2):box(4)) = 0;
imgRGB(box(1):box(3),box(2):box(4),1) = 0;
imgRGB(box(1):box(3),box(2):box(4),2) = 0;
imgRGB(box(1):box(3),box(2):box(4),3) = 0;

% get edges
edges = uint8(edgesBefore*255);
edges(box(1):box(3),box(2):box(4)) = 100;

% get inducers.
% merge adjacent edge points (assuming there aren't more than 2 adjacent edge pixels)
inducers = [];
y=box(1)-1;
for x=box(2)-1:box(4)+1
    % if not edge
    if edges(y,x) ~= 255
        continue;
    end
    
    if edges(y,x) ~= edges(y,x+1)
        % if edge pixel has no adjacent edge pixels take it as is
        inducers = [inducers; [x,y]];
    else
        % if edge pixel is followed by another one, take their middle
        inducers = [inducers; [x+0.5,y]];
        edges(y,x+1) = 0;
    end
end
y=box(3)+1;
for x=box(2)-1:box(4)+1
    % if not edge
    if edges(y,x) ~= 255
        continue;
    end
    
    if edges(y,x) ~= edges(y,x+1)
        % if edge pixel has no adjacent edge pixels take it as is
        inducers = [inducers; [x,y]];
    else
        % if edge pixel is followed by another one, take their middle
        inducers = [inducers; [x+0.5,y]];
        edges(y,x+1) = 0;
    end
end
x=box(2)-1;
for y=box(1)-1:box(3)+1
    % if not edge
    if edges(y,x) ~= 255
        continue;
    end
    
    if edges(y,x) ~= edges(y+1,x)
        % if edge pixel has no adjacent edge pixels take it as is
        inducers = [inducers; [x,y]];
    else
        % if edge pixel is followed by another one, take their middle
        inducers = [inducers; [x,y+0.5]];
        edges(y+1,x) = 0;
    end
end
x=box(4)+1;
for y=box(1)-1:box(3)+1
    % if not edge
    if edges(y,x) ~= 255
        continue;
    end
    
    if edges(y,x) ~= edges(y+1,x)
        % if edge pixel has no adjacent edge pixels take it as is
        inducers = [inducers; [x,y]];
    else
        % if edge pixel is followed by another one, take their middle
        inducers = [inducers; [x,y+0.5]];
        edges(y+1,x) = 0;
    end
end
inducers = double(inducers);
if size(inducers,1) < 2
    disp(['Error - only ' num2str(size(inducers,1)) ' inducers'])
end

% get inducer orientations.
for i=1:size(inducers,1)
    p = inducers(i,:);
    % get gradient vector, normalize it, flip y axis
    gd = [edgeOrDX(p(2),p(1)),edgeOrDY(p(2),p(1))];
    gd = gd/norm(gd);
    gd(2) = -gd(2);
    % get angle between the gradient vector and the x axis
    or = atan2(gd(2),gd(1));
    if or<0
        or = or + 2*pi;
    end
    % convert to edge direction
    or = mod(or + pi/2, 2*pi);
    
    inducerOrs(i) = or;
end

% direct orientation to the box. an edge to the left should have
% orientation between -pi/2 and pi/2, etc.
for i=1:size(inducers,1)
    x = inducers(i,1);
    y = inducers(i,2);
    or = inducerOrs(i);
    if x < box(2) % left
        if or > pi/2 && or < 3*pi/2
            or = mod(or + pi,2*pi);
        end
    elseif x > box(4) % right
        if or < pi/2 || or > 3*pi/2
            or = mod(or + pi,2*pi);
        end
    elseif y < box(1) % above
        if or < pi
            or = mod(or + pi,2*pi);
        end
    else % below
        if or > pi
            or = mod(or + pi,2*pi);
        end
    end
    inducerOrs(i) = or;
end

% flip y axis (so y axis points up)
inducers(:,2) = size(img,2) - inducers(:,2) + 1;

% complete curves.
% for each pair complete a curve between them and calculate the probability
% for that curve. connect inducer curves greedily according to the most
% probable curves.

indPairs = nchoosek(1:size(inducers,1),2);
indPairProb = zeros(size(indPairs,1),1);
indPairCurve = cell(size(indPairs,1),1);

for i=1:size(indPairs,1)
    ind1 = indPairs(i,1);
    ind2 = indPairs(i,2);
    
    relEndPoint = transPoints(inducers(ind2,:), inducers(ind1,:), inducerOrs(ind1));
    % ------
    % get mean curve
    
    endPoint = relEndPoint;
    endPointBin = floor((endPoint - [relMinX, relMinY])/binSize) + 1;
    
    % get curve probability
    flippedY = size(numSamples,2) - endPointBin(2) + 1;
    curveProb = numSamples(flippedY, endPointBin(1));
    
    allEndPoints = [[frags{:,1}]', [frags{:,2}]'];
    endPointFrags = find(allEndPoints(:,1)==endPointBin(1) & allEndPoints(:,2)==endPointBin(2));
    
    numFrags = numel(endPointFrags);
    if numFrags < 1
%         disp(['No sample curves for pair' num2str(i)]);
        continue;
    end
    
    allRepPts = zeros(numFrags,numCurveRepPts*2); % all curves' representative points
    % times 2 because each point is x,y
    for j=1:numFrags
        fragPts = frags{endPointFrags(j),3};
        fragPts = [[0,0]; fragPts];
        
        repPts = getCurveEquiPoints(fragPts, numCurveRepPts);
        allRepPts(j,:) = reshape(repPts,1,numCurveRepPts*2);
    end
    
    % calculate mean curve
    %         coeff = pca(allRepPts);
    meanPts = mean(allRepPts,1);
    meanPts = reshape(meanPts, numCurveRepPts, 2);
    meanCurve = [[0,0]; meanPts; endPoint];
    % scatter(meanPts(:,1),meanPts(:,2),7,'r','filled');
    % line(meanPts(:,1),meanPts(:,2),'Color','r');
    
    % ------
    
    % get completion curve (image coordinates)
    completion = meanCurve;
    % rotate back to the starting point's original orientation and return to image location
    completion = transBackPoints(completion, inducers(ind1,:), inducerOrs(ind1));
    % return to image location
%     completion = completion + repmat(inducers(ind1,:),size(completion,1),1);
    % flip y axis (so y axis points down)
    completion(:,2) = size(img,2) - completion(:,2) + 1;
    
    indPairCurve{i} = completion;
    indPairProb(i) = curveProb;
end

% choose curves
newCurves = cell(0,1);
numNewCurves = 0;
while any(indPairProb>0)
    [~,i]=max(indPairProb);
    ind1 = indPairs(i,1);
    ind2 = indPairs(i,2);
    
    numNewCurves = numNewCurves + 1;
    newCurves{numNewCurves,1} = indPairCurve{i};
    
    indPairProb(indPairs(:,1)==ind1) = 0;
    indPairProb(indPairs(:,1)==ind2) = 0;
    indPairProb(indPairs(:,2)==ind1) = 0;
    indPairProb(indPairs(:,2)==ind2) = 0;
end

% imshow(imgRGB)
% hold on
% line(c1(:,1),c1(:,2),'Color','r');
% figure
edges(edgesBefore) = 255;
imshow(edgesBefore)
figure
imshow(imgRGB)
set(gcf,'position',get(0,'screensize'))

% plot curves
for i=1:size(newCurves,1)
    c = newCurves{i};
    
    hold on
    line(c(:,1),c(:,2),'Color','r');
    hold on
    scatter(c(1,1),c(1,2),7,'y','filled');
end
