%% Generate test set curves for benchmark

numNeededSamples = 5000; % number of sample curves to collect

% when randomly choosing the curve part to leave occluded, keep some unused
% points on each side to allow the calculation of orientation at the end
% points
numMarginPts = 5;

% do not use curves with too little points
shortestCurvePts = 5;

% do not use curves that are too short (Eucledean distance)
shortestCurveEuc = 5;

% calculate orientation by fitting a line using points at most gapSize away
gapSize = 3;

% place curves in bins according to their size to make sure with have a
% uniform distribution of curves by size
binsMaxSize = 200;
numBins = 10;
binSize = binsMaxSize / numBins;
numCurvesPerBin = floor(numNeededSamples / numBins);

% set folder paths
params.datasetFolder = '../data/curve fragments dataset/CFGD_release/';
params.curvesFolder = [params.datasetFolder 'GT_mat_CFGD_format/'];
params.imgsFolder = [params.datasetFolder 'img/'];
params.outFolder = '../data/';

% test images subset percent
testSetPercent = 0.1;

% get all dataset image names and sizes
files = dir([params.imgsFolder '*jpg']);
numImgs = length(files);
imgNames = cell(size(files,1),1);
imgSizes = zeros(size(files,1),2);
for i=1:size(files,1)
    imgNames{i} = files(i).name;
    img = imread([params.imgsFolder files(i).name]);
    imgSizes(i,:) = [size(img,1), size(img,2)];
end

% get random test images 
numTestImgs = floor(numImgs * testSetPercent);
testImgs = randperm(numImgs, numTestImgs);
testImgNames = {};
for i = 1:length(testImgs)
    testImgNames{i} = imgNames{testImgs(i)};
end

% count total number of curves and curve lengts (unneeded)
numCurves = 0;
curveLens = [];
for i = testImgs
    % go over the different annotators
    for annotatorNum=1:3
        % load image curves
        imgName = imgNames{i};
        baseName = imgName(1:end-4);
        data = load([params.curvesFolder baseName],'groundTruth');
        curves = data.groundTruth{annotatorNum}; % use one of the annotators

        numCurves = numCurves + length(curves);
        for j=1:length(curves)
%             curveLen = size(curves{j},1);
            curveLen = sqrt(sum((curves{j}(1,1:2) - (curves{j}(end,1:2))) .^ 2));
            curveLens = [curveLens; curveLen];
        end
    end
end
disp(['num curves = ' num2str(numCurves)]);
% hist(curveLens)


% collect curves
testCurves = {};
testCurvesMat = zeros(6,0);
numBinCurves = zeros(numBins, 1); % count curves in each bin
numCurves = 0; % count total number of curves
isDone = false;
while ~isDone
    [min(numBinCurves), mean(numBinCurves)]
    
    % randomly pick image
    i = randi(length(testImgs));
    imgNum = testImgs(i);
    % randomly pick annotator
    annotatorNum = randi(3);
    
    % load image curves
    imgName = imgNames{imgNum};
    baseName = imgName(1:end-4);
    imgSize = imgSizes(imgNum,:);
    data = load([params.curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{annotatorNum}; % get curves (0-based coords)
    
    % randomly pick curve
    curveNum = randi(length(curves));
    c = curves{curveNum};
    c = fixCurve(c, imgSize);
    
    % remove points outside image boundary
    outPts = logical(c(:,1)<1 | c(:,2)<1 | c(:,1)>imgSizes(imgNum,2) | c(:,2)>imgSizes(imgNum,1));
    c = c(~outPts, :);
    
    % randomly flip curve direction
    isFlipped = 0;
    if rand > 0.5
        c = c(size(c,1):-1:1,:);
        isFlipped = 1;
    end
    
    ptsWithMargin = c(:,1:2);
    if size(ptsWithMargin, 1) <= 2*numMarginPts + shortestCurvePts
        continue;
    end
    pts = ptsWithMargin(numMarginPts+1:size(ptsWithMargin, 1)-numMarginPts, :);
    numPts = size(pts, 1);
    
    % pick number of points to use as occluded part
    numOccPts = randi(numPts-shortestCurvePts) + shortestCurvePts;
    
    % pick start point and infer end point
    p1 = randi(numPts - numOccPts + 1);
    p2 = p1 + numOccPts - 1;
    
    % update appropriate bin
    curveLen = sqrt(sum((pts(p1,:) - pts(p2,:)) .^ 2));
    curveLenBin = floor(curveLen / binSize) + 1;
    if curveLenBin > numBins
        curveLenBin = numBins;
    end
    
    % ignore too short curves (in number of pixels)
    if curveLen < shortestCurveEuc
        continue;
    end
    
    % don't collect more than needed
    if numBinCurves(curveLenBin) >= numCurvesPerBin
        continue;
    end
    
    % add margin
    p1 = p1 + numMarginPts;
    p2 = p2 + numMarginPts;
    
    % calculate inducer orientations - if available use margin points such
    % that the arc-length for the calculation of orienation would be equal
    % to gapSize
    arcLenToP1 = getArcLength(ptsWithMargin,p1,1);
    arcLenFromP2 = getArcLength(ptsWithMargin,p2,size(ptsWithMargin,1));
    p1o = find(arcLenToP1 < gapSize,1);
    p2o = find(arcLenFromP2 < gapSize,1,'last');
    if p1o == p1
        p1o = p1 - 1;
    end
    if p2o == p2
        p2o = p2 + 1;
    end
    p1Or = getOrFit(ptsWithMargin(p1o:p1,:));
    p2Or = getOrFit(ptsWithMargin(p2o:-1:p2,:));
    
    % add to struct
    numCurves = numCurves + 1;
    p1xy = ptsWithMargin(p1,:);
    p2xy = ptsWithMargin(p2,:);
    numBinCurves(curveLenBin) = numBinCurves(curveLenBin) + 1;
    testCurves{numCurves}.pts = ptsWithMargin; % all curve points
    testCurves{numCurves}.p1 = p1; % index of first inducer
    testCurves{numCurves}.p2 = p2; % index of second inducer
    testCurves{numCurves}.p1Or = p1Or; % orientation of first inducer
    testCurves{numCurves}.p2Or = p2Or; % orientation of second inducer
    testCurves{numCurves}.imgName = imgName;
    testCurves{numCurves}.annotatorNum = annotatorNum;
    testCurves{numCurves}.curveNum = curveNum;
    testCurves{numCurves}.isFlipped = isFlipped; % if curve was flipped
    
    % add to mat
    testCurvesMat = [testCurvesMat; [p1xy(1), p1xy(2), p1Or, p2xy(1), p2xy(2), p2Or]];
    
    if min(numBinCurves) >= numCurvesPerBin
        isDone = true;
    end
end

testSet.curves = testCurves;
testSet.imgNames = testImgNames;
testSet.imgNums = testImgs;

save([params.outFolder 'test_set/testSet.mat'], 'testSet');
dlmwrite([params.outFolder 'test_set/testCurvesMat.txt'], testCurvesMat);

disp('done')


% Show histogram of lengths

lens = zeros(length(testCurves), 2);
for i = 1:length(testCurves)
    lens(i,1) = testCurves{i}.p2 - testCurves{i}.p1 + 1;
    lenEuc = sqrt(sum((testCurves{i}.pts(testCurves{i}.p1,:) - testCurves{i}.pts(testCurves{i}.p2,:)) .^ 2));
    lens(i,2) = lenEuc;
end

hist(lens(:,1))
title('Histogram of curve lengths (number of points)');
figure
hist(lens(:,2))
title('Histogram of curve lengths (distance between start and end points)');


%% Visualize curves


for i = 1:numel(testCurves)
    % curve points
    pts = testCurves{i}.pts(testCurves{i}.p1:testCurves{i}.p2,:);
    % inducers
    p1xy = testCurves{i}.pts(testCurves{i}.p1,:);
    p2xy = testCurves{i}.pts(testCurves{i}.p2,:);
    p1Or = testCurves{i}.p1Or;
    p2Or = testCurves{i}.p2Or;
    
    c=pts;
    line(c(:,1),c(:,2),'LineWidth',2);
    hold on
    visInducers(p1xy, p1Or, p2xy, p2Or);
    axis([1 481 1 481]);
    saveas(gcf,[params.outFolder 'test_set/curve_imgs/' num2str(i) '.png']);
    close all
end

