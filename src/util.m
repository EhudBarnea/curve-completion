
%% params


if ~exist('completeCurve.m','file')
    disp('Error: file completeCurve.m cannot be found, cd to the src folder');
end

addpath('kde2d/');

% set folder paths
params.datasetFolder = '../data/curve fragments dataset/CFGD_release/';
params.curvesFolder = [params.datasetFolder 'GT_mat_CFGD_format/'];
params.imgsFolder = [params.datasetFolder 'img/'];
params.outFolder = '../data/';

% number of parpool workers
param.numParWorkers = 1; % default
% param.numParWorkers = 20; % sge44
% param.numParWorkers = 32; % sge59

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
    params.matchDistFactor = 10;
else
    params.matchDistFactor = 1;
end
params.matchOr = 10*pi/180;
% params.matchOr = 2*pi/180;
params.matchSI = true; % match in a scale invariant way

% set jet colormap
% set(groot,'DefaultFigureColormap',jet)

%% train / collect data

collectCurveFrags(params);
load([params.outFolder 'all_frags/frags']);
collectScaleInvCurveFrags(frags, params);

%% complete a single curve

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags']);
end
if ~exist('fragsSI','var') && params.matchSI
    load([params.outFolder 'all_frags/fragsSI']);
end
 
vis = true;
%--
% p1 = [30,30];
% or1 = 0;
% p2 = [30,90];
% or2 = 0;
%--
p1 = [-30,0];
or1 = deg2rad(30);
p2 = [30,0];
or2 = pi-or1;
%--


if ~params.matchSI
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, vis);
else
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, fragsSI, params, vis);
end


%% run completion demo

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags']);
end
if ~exist('fragsSI','var') && params.matchSI
    load([params.outFolder 'all_frags/fragsSI']);
end
 
%img = imread('a.png');
demoCompletionSingle(frags, img, params);

%% show completions of all curves

if ~exist('frags','var') && ~params.matchSI
    load([params.outFolder 'all_frags/frags']);
end
if ~exist('fragsSI','var') && params.matchSI
    load([params.outFolder 'all_frags/fragsSI']);
end
 

if ~params.matchSI
    completeAllCurves(frags, params);
else
    completeAllCurves(fragsSI, params);
end

%% visualize number of samples (non scale invariant)

if params.matchSI
    disp('Error: incorrect configuration');
end
if ~exist('frags','var')
    load([params.outFolder 'all_frags/frags']);
end
for i=1:8
    
    % flip numSamples to image [y,x] coordinates
    tmp = numSamples(:,:,i);
    tmp = tmp';
    tmp = flipud(tmp);
    
    imagesc(tmp>20);
    saveas(gcf,[params.outFolder '/tmp/numSamples' num2str(i) '.png']);
    close all
end

%% Analyse scale invariance - single configuration

if params.matchSI
    disp('Error: incorrect configuration');
end
if ~exist('frags','var')
    load([params.outFolder 'all_frags/frags']);
end

visEachScale = true;

endPointDirection = deg2rad(340);
endPointOr = deg2rad(135);
analyzeScaleInvariance(endPointDirection, endPointOr, frags, params, visEachScale);
display('done');

%% Analyse scale invariance - all configurations

if params.matchSI
    disp('Error: incorrect configuration');
end
if ~exist('frags','var')
    load([params.outFolder 'all_frags/frags']);
end

visEachScale = true;

% loop over different end point directions and end point orientations. The range variables are used since
% parfor requires consecutive integers.
rangeDir = deg2rad(270:5:360);
rangeOr = deg2rad(0:45:315);
% get all pairs of endPointDirection and endPointOr
[p,q] = meshgrid(rangeDir, rangeOr);
rangePairs = [p(:) q(:)];
parpool(param.numParWorkers);
parfor i = 1:size(rangePairs,1)
    disp([num2str(i) '/' num2str(size(rangePairs,1)) ' start']);
    endPointDirection = rangePairs(i,1);
    endPointOr = rangePairs(i,2);
    analyzeScaleInvariance(endPointDirection, endPointOr, frags, params, visEachScale);
    disp([num2str(i) '/' num2str(size(rangePairs,1)) ' done']);
end

