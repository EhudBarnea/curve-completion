
%% params


if ~exist('completeCurve.m','file')
    disp('Error: file completeCurve.m cannot be found, cd to the src folder');
end

% addpath('kde2d/');
addpath('DiscreteFrechetDist');

% set folder paths
params.datasetFolder = '../data/curve fragments dataset/CFGD_release/';
params.curvesFolder = [params.datasetFolder 'GT_mat_CFGD_format/'];
params.imgsFolder = [params.datasetFolder 'img/'];
params.outFolder = '../data/';

% test set file
params.testSetFilename = [params.outFolder 'test_set/testSet.mat'];

% get all dataset image names and sizes
files = dir([params.imgsFolder '*jpg']);
params.numImgs = length(files);
params.imgNames = cell(size(files,1),1);
params.imgSizes = zeros(size(files,1),2);
for i=1:size(files,1)
    params.imgNames{i} = files(i).name;
    img = imread([params.imgsFolder files(i).name]);
    params.imgSizes(i,:) = [size(img,1), size(img,2)];
end

% ----- params for the "frags" dataset
% all the curve fragments in cannonical pose are kept in spatio-angular
% bins according to these parameters

% minimum distance between two points for which we keep the curve
% fragments.
params.minDist = 0; % in pixels (Eucledean)

params.binSize = 1; % in pixels
% size relative to reference point (adjusted to be a multiple of binSize)
params.relMaxX = 200;
params.relMaxY = 200;
params.relMaxX = floor(params.relMaxX/params.binSize)*params.binSize;
params.relMaxY = floor(params.relMaxY/params.binSize)*params.binSize;
params.relMinX = -params.relMaxX;
params.relMinY = -params.relMaxY;
params.numBins = [(params.relMaxX-params.relMinX)/params.binSize, (params.relMaxY-params.relMinY)/params.binSize];

% discretization of end point orientation
params.numOrBins = 8;
params.orBinSize = 2*pi/params.numOrBins;
% -----

% number of points in the generated curve
params.numCurveRepPts = 30;

% ----- params for the "fragsSI" dataset
params.siNumAngBins = 30; % Number of angular bins for curve end points
params.siAngBinSize = pi/params.siNumAngBins;
params.siNumOrBins = 8;
params.siOrBinSize = 2*pi/params.siNumOrBins;
% -----
params.siMinScale = 20; % do not use curve of scale smaller than this

% match distance factor and orientation. This configures which curve fragments to
% use (the close enough fragments) for a completion given start and end
% inducers in cannonical pose. The distance factor is in number of pixels,
% but relative to the distance between the start and end inducers of a
% curve fragment (if relMatchDist=true).
% For example, a matchDistFactor of 10 means that for curves for which the
% (Euclidean) distance between the inducers is 10 pixels, other curve
% fragments will be considered if their end point is 1 pixel away.
% set to a negative value to disable (using only curves that end in the
% exact same pixel
params.relMatchDist = true;
if params.relMatchDist
    params.matchDistFactor = 8;
else
    params.matchDistFactor = 1;
end
params.matchOr = 2*pi/180;

% maximum number of fragments to use (less = faster)
params.maxFragsToUse = inf;

params.matchSI = true; % match in a scale invariant way

% params for curve generation assuming extensibility
params.extensible = true;
params.numMinFrags = 400;
% params.numMinFrags = inf;
params.extMatchOr = params.matchOr;

%% train / collect data

% create structures (containing both train and test data)
collectCurveFrags(params);
load([params.outFolder 'all_frags/frags.mat']);
collectScaleInvCurveFrags(frags, [], params);

% create structure for training data
testImgNums = [2, 10, 7, 28, 27];
trainImgNums = setdiff(1:params.numImgs, testImgNums);
collectScaleInvCurveFrags(frags, trainImgNums, params);

%% complete a single curve

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags.mat']);
end
if ~exist('fragsSI','var') && params.matchSI
    load([params.outFolder 'all_frags/fragsSI.mat']);
end

vis = true;

opt = 2;
if opt == 1
    p1 = [30,30];
    or1 = 0;
    p2 = [30,90];
    or2 = 0;
elseif opt == 2
    p1 = [-30,0];
    or1 = deg2rad(30);
    p2 = [30,0];
    or2 = pi-or1;
elseif opt == 3
    p1 = [0,0];
    or1 = 0;
    p2 = [40,-80];
    or2 = pi/2;
end


if ~params.matchSI
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, vis);
else
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, fragsSI, params, vis);
end


%% run completion demo

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags.mat']);
end
if ~exist('fragsSI','var') && params.matchSI
    load([params.outFolder 'all_frags/fragsSI.mat']);
end

% filename = '../data/test_imgs/https-::commons.wikimedia.org:wiki:File-Mushroom-IMG_3304.JPG';
% filename = '../data/test_imgs/hipo.jpg';
% filename = '../data/test_imgs/http:::pngimg.com:img:flowers:camomile.png';

filename = '../data/test_imgs/screw2.jpg';
% filename = '../data/test_imgs/flower1.png';

img = imread(filename);
if ~params.matchSI
    demoCompletionSingle(frags, img, params);
else
    demoCompletionSingle(fragsSI, img, params);
end

%% show completions of all curves

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags.mat']);
end
if ~exist('fragsSI','var') && params.matchSI
    load([params.outFolder 'all_frags/fragsSI.mat']);
end


if ~params.matchSI
    resCurves = completeAllCurves(frags, params);
else
    resCurves = completeAllCurves(fragsSI, params);
end

%% run completion over test set

% params.maxFragsToUse = 1;
params.maxFragsToUse = 600;
% params.extensible = false;
params.extensible = true;
vis_res = false;

% Load train curve fragments
if ~params.matchSI
    disp('Error: not implemented');
end
if ~exist('fragsSI','var')
    load([params.outFolder 'all_frags/fragsSI_train.mat']);
end

% load test set
data = load(params.testSetFilename,'testSet');
testCurves = data.testSet.curves;

% complete test curves curves
testCompletions = cell(size(testCurves));
for i = 1:numel(testCurves)
    i
    %     tic
    
    % get inducers
    p1 = testCurves{i}.pts(testCurves{i}.p1, :);
    p2 = testCurves{i}.pts(testCurves{i}.p2, :);
    or1 = testCurves{i}.p1Or;
    or2 = testCurves{i}.p2Or;
    
    % complete curve
    [c, ~, ~] = completeCurve(p1, or1, p2, or2, fragsSI, params, false);
    
    % visualize ground-truth vs completion
    if vis_res
        % get gt curve
        gt = testCurves{i}.pts(testCurves{i}.p1:testCurves{i}.p2,:);
        
        % calculate reconstruction error
        [error, relError] = curve_dist(gt, c);
        
        line(c(:,1),c(:,2),'color','blue','LineWidth',2);
        hold on
        line(gt(:,1),gt(:,2),'color','k','LineWidth',2);
        hold on
        visInducers(p1, or1, p2, or2);
        axis([1 481 1 481]);
        title(['Error = ' num2str(error) ' Relative error = ' num2str(relError)])
        saveas(gcf,[params.outFolder 'test_set/gt_vs_completion/' num2str(i) '.png']);
        close all
    end
    
    testCompletions{i} = c;
    %     toc
end

% save completions
save([params.outFolder 'test_set/testSet_completions_meancurve.mat'], 'testCompletions');
% save([params.outFolder 'test_set/testSet_diff_completions_meancurve.mat'], 'testCompletions');

%% Evaluate completions as test set reconstructions

% load test set
data = load(params.testSetFilename,'testSet');
testCurves = data.testSet.curves;

% evaluate spiral reconstruction results
numSpiralPts = 200;
testCompletions_spiralMat = dlmread([params.outFolder 'test_set/testSet_completions_spiral.txt']);
% testCompletions_spiralMat = dlmread([params.outFolder 'test_set/testSet_diff_completions_spiral.txt']);
testCompletions_spiral = cell(size(testCurves));
for i=1:length(testCurves)
    s = 2+(1+numSpiralPts)*(i-1);
    e = s+200-1;
    testCompletions_spiral{i} = testCompletions_spiralMat(s:e,:);
end
[spiral_aarc, spiral_errors, spiral_sortID] = evalReconst(testCompletions_spiral, testCurves, 'b');
hold on
%
% evaluate mean curve reconstruction results
data = load([params.outFolder 'test_set/testSet_completions_meancurve.mat'], 'testCompletions');
% data = load([params.outFolder 'test_set/testSet_diff_completions_meancurve.mat'], 'testCompletions');
testCompletions = data.testCompletions;
[mc_aarc, mc_errors, mc_sortID] = evalReconst(testCompletions, testCurves, 'r');


% add legend
lgd = legend('Euler spiral', 'Mean curve');
lgd.FontSize = 16;

[spiral_aarc, mc_aarc]

disp('done')

%% Visualize completions over test set

lineWidth = 1;

% load test set
data = load(params.testSetFilename,'testSet');
testCurves = data.testSet.curves;

% load Euler spiral reconstruction results
numSpiralPts = 200;
testCompletions_spiralMat = dlmread([params.outFolder 'test_set/testSet_completions_spiral.txt']);
testCompletions_spiral = cell(size(testCurves));
for i=1:length(testCurves)
    s = 2+(1+numSpiralPts)*(i-1);
    e = s+200-1;
    testCompletions_spiral{i} = testCompletions_spiralMat(s:e,:);
end

% load mean curve reconstruction results
data = load([params.outFolder 'test_set/testSet_completions_meancurve.mat'], 'testCompletions');
% data = load([params.outFolder 'test_set/testSet_diff_completions_meancurve.mat'], 'testCompletions');
testCompletions_meanCurve = data.testCompletions;

for i = 1:numel(testCurves)
    i
    
    % get image
    imgName = [params.datasetFolder 'img/' testCurves{i}.imgName];
    img = imread(imgName);
    imshow(img)
    
    % get gt curve
    gt = testCurves{i}.pts(testCurves{i}.p1:testCurves{i}.p2,:);
    % flip y axis
    gt(:,2) = size(img,1) - gt(:,2) + 1;
    
    % get inducers
    p1 = testCurves{i}.pts(testCurves{i}.p1, :);
    p2 = testCurves{i}.pts(testCurves{i}.p2, :);
    or1 = testCurves{i}.p1Or;
    or2 = testCurves{i}.p2Or;
    % flip y axis
    p1(2) = size(img,1) - p1(2) + 1;
    p2(2) = size(img,1) - p2(2) + 1;
    scale = int32(norm(p1 - p2));
    
    % get completions
    c_spiral = testCompletions_spiral{i};
    c_meanCurve = testCompletions_meanCurve{i};
    % flip y axis
    c_spiral(:,2) = size(img,1) - c_spiral(:,2) + 1;
    c_meanCurve(:,2) = size(img,1) - c_meanCurve(:,2) + 1;
    
    % calculate reconstruction error
%     [error, relError] = curve_dist(gt, c_spiral);
%     title(['Error = ' num2str(error) ' Relative error = ' num2str(relError)])
    
    line(gt(:,1),gt(:,2),'color','magenta','LineWidth',lineWidth);
    hold on
    line(c_spiral(:,1),c_spiral(:,2),'color','blue','LineWidth',lineWidth);
    hold on
    line(c_meanCurve(:,1),c_meanCurve(:,2),'color','green','LineWidth',lineWidth);
    hold on
    
%     visInducers(p1, or1, p2, or2);
    saveas(gcf,[params.outFolder 'test_set/gt_vs_completions/' num2str(scale) '_' num2str(i) '.png']);
    close all
end


%% visualize number of samples (non scale invariant)

if params.matchSI
    disp('Error: incorrect configuration');
end
if ~exist('frags','var')
    load([params.outFolder 'all_frags/frags.mat']);
end
for i=1:8
    
    % flip numSamples to image [y,x] coordinates
    tmp = numSamples(:,:,i);
    tmp = tmp';
    tmp = flipud(tmp);
    tmp = tmp(1+size(tmp,1)/2:end,:);
    
    tmp(tmp>300) = 300;
    imagesc(-tmp);
    %     imagesc(log(tmp));
    %     imagesc(tmp>200);
    colormap gray
    hcb=colorbar;
    tcks = -hcb.Ticks;
    
    tcks=num2cell(tcks);
    tcks{1}=['>' num2str(tcks{1})];
    set(hcb,'TickLabels',tcks)
    set( hcb, 'YDir', 'reverse' );
    axis equal
    xlim([0,400])
    ylim([0,300])
    saveas(gcf,[params.outFolder '/tmp/numSamples' num2str(i) '.png']);
    close all
end

%% Analyse scale invariance - all configurations

params.matchSI=false;
params.extensible=false;

if ~exist('frags','var')
    load([params.outFolder 'all_frags/frags.mat']);
end

visEachScale = false;

measuresOutFolder = [params.outFolder 'stats/scale_analysis/'];



% loop over different end point orientations
pr1Ors = deg2rad(0:22.5:180);
pr2Ors = deg2rad(0:22.5:337.5);
% get all pairs
[p,q] = meshgrid(pr2Ors, pr1Ors);
% prepare matrices to keep analysis results
confMeanVar = zeros(size(p));
confMeanDiffMu = zeros(size(p));
% confMeanDiffSTD = zeros(size(p));
orPairs = [p(:) q(:)];
% for i = 1:size(orPairs,1)
% parpool(20)
parfor i = 1:size(orPairs,1)
    fprintf('%d/%d start\n',i,size(orPairs,1));
    
    % get p2 and or2 in relative to p1 and or1
    [endPoint, endPointOr] = transPoints([1,0], orPairs(i,1), [0,0], orPairs(i,2));
    
    
    % analyze scales
    [meanVar, meanDiffMu, meanDiffVar] = analyzeScaleInvariance(endPoint, endPointOr, frags, params, visEachScale);
    confMeanVar(i) = meanVar;
    confMeanDiffMu(i) = meanDiffMu;
    %     confMeanDiffVar(i) = meanDiffVar;
    
    fprintf('%d/%d done\n',i,size(orPairs,1));
end
usedParams = params;
save([measuresOutFolder 'measures.mat'],'confMeanVar','confMeanDiffMu','usedParams');

%% Visualize analysis above
imagesc(confMeanDiffMu/100)
colorbar
colormap jet
set(gca,'FontSize',22,'linewidth',2)
yticks([1, 3, 5, 7, 9]);
yticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi',});
% xticks([1, 3, 5, 7, 9, 11, 13, 15]);
% xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4'});
xticks([1, 5, 9, 13, 15]);
xticklabels({'0','\pi/2','\pi','3\pi/2'});
set(gcf, 'Color', 'w');
xlhand = get(gca,'xlabel');
set(xlhand,'string','\theta_2','fontsize',37);
ylhand = get(gca,'ylabel');
set(ylhand,'string','\theta_1','fontsize',37,'rot',0,'Position',get(ylhand,'Position') - [0.3 0 0]);
% set(gca,'FontSize',13,'FontWeight','bold','linewidth',2)
scale = 0.03;
pos = get(gca, 'Position');
pos(1) = pos(1)+scale*pos(3);
pos(3) = (1-scale)*pos(3);
set(gca, 'Position', pos)
% set(gcf, 'Units', 'pixels', 'Position', [10, 100, 700, 400]);
export_fig([params.outFolder 'stats/scale_analysis/meanCurveVariance.pdf']);
close all
%
imagesc(confMeanVar/100)
colorbar
colormap jet
set(gca,'FontSize',22,'linewidth',2)
yticks([1, 3, 5, 7, 9]);
yticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'});
% xticks([1, 3, 5, 7, 9, 11, 13, 15]);
% xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4'});
xticks([1, 5, 9, 13, 15]);
xticklabels({'0','\pi/2','\pi','3\pi/2'});
set(gcf, 'Color', 'w');
xlhand = get(gca,'xlabel');
set(xlhand,'string','\theta_2','fontsize',37);
ylhand = get(gca,'ylabel');
set(ylhand,'string','\theta_1','fontsize',37,'rot',0,'Position',get(ylhand,'Position') - [0.3 0 0]);
scale = 0.03;
pos = get(gca, 'Position');
pos(1) = pos(1)+scale*pos(3);
pos(3) = (1-scale)*pos(3);
set(gca, 'Position', pos)
% set(gcf, 'Units', 'pixels', 'Position', [10, 100, 700, 400]);
export_fig([params.outFolder 'stats/scale_analysis/fragVariances.pdf']);
close all



%% Generate summary figure


params.matchSI = false;
params.extensible = false;

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags.mat']);
end


vis = true;

s = 1.7;
p1 = [-30,0]*s;
or1 = deg2rad(30);
p2 = [30,0]*s;
or2 = pi-or1;


[c, isUsable, out] = completeCurve_paperFig(p1, or1, p2, or2, frags, params, vis);

%% Generate reconstruction figure
% to run, disable 'figure' and 'imshow' commands in 'demoCompletionSingle.m'

% params.extensible = true;
if ~exist('fragsSI','var')
    load([params.outFolder 'all_frags/fragsSI.mat']);
end

opt =2;
c=vec2;
xx=0;
if opt==1
    filename = '../data/test_imgs/screw2.jpg';
    xx=30;
    inducers = [151.1983  187.7905    2.9000  303.6844  488.6285    2.5400];
end
if opt==2
    filename = '../data/test_imgs/flower.png';
    %     inducers = [470.0785  217.4611    5.0854  583.3741   36.1881    3.9200];
    inducers = [587.9059   39.5870    4.0921  467.8126  223.1259    5.0612];
end
if opt==3
    filename = '../data/test_imgs/mushroom.jpg';
    inducers = [459.3    317.3    5.7    1661.8    406.3    3.85];
end
if opt==4
    filename = '../data/test_imgs/mushroom2.jpg';
    inducers = [293.6927  301.4071    5.4146  880.3352  267.6920    3.8967];
end
if opt==5
    filename='../data/test_imgs/smiley.jpg';
    inducers = [157.0378  145.0107    5.8650  108.1891  174.9860    0.2111];
end

img = imread(filename);
img = imtranslate(img,[xx, 0],'FillValues',255);
% img(100:406,459:1661,:) = 255;
imshow(img)
line(c(:,1),c(:,2),'Color','magenta', 'LineWidth',6)

hold on
% demoCompletionSingle(fragsSI, img, params);
demoCompletionSingle(fragsSI, img, params, inducers);
saveas(gcf,['/Users/ehud/Downloads/tmp/' num2str(opt) '.png']);
close all
%%
for i=1%:2000
    i
    inducers2=inducers;
    %     inducers2(6) = inducers(6)-i*0.01;
    
    img = imread(filename);
    imshow(img)
    line(c(:,1),c(:,2),'Color','magenta', 'LineWidth',3)
    hold on
    demoCompletionSingle(fragsSI, img, params, inducers2);
    saveas(gcf,['/Users/ehud/Downloads/tmp/' num2str(i) '.png']);
    close all
end

%% Analyze extent of extensibility
% compare all completions with scale invariance and with scale invariance
% but without extensibility (when enough fragments were available). Report
% the maximal Frechet distance between curves.

params.extensible=false;
resCurves_si = completeAllCurves(fragsSI, params);
params.extensible=true;
params.numMinFrags = inf;
resCurves_si_ext = completeAllCurves(fragsSI, params);

%%
minFrags = 4000;

errors = [];
for minFrags = 0:1:3000
    minFrags
    res = [];
    for i = 1:size(resCurves_si, 2)
        c1 = resCurves_si{i}.pts;
        c2 = resCurves_si_ext{i}.pts;
        n1 = resCurves_si{i}.numFrags;
        n2 = resCurves_si_ext{i}.numFrags;
        
        if n1 < minFrags || n2 < minFrags
            continue;
        end
        
        [d, dRel] = curve_dist(c1, c2);
        res = [res; [d, dRel]];
    end
    % res
%     max(res)
%     mean(res)
%     median(res)
    errors = [errors; [minFrags, max(res), mean(res)]];
end

% results for minFrags=400:
% max =    9.8471    0.0492
% mean =   5.0291    0.0251
% median = 4.4428    0.0222

% results for minFrags=2000:
% max =    3.5873    0.0179
% mean =   3.4325    0.0172
% median = 3.4203    0.0171


