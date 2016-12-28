
%% params


addpath('export_fig/');

params.datasetFolder = '../data/curve fragments dataset/CFGD_release/';
params.curvesFolder = [params.datasetFolder 'GT_mat_CFGD_format/'];
params.imgsFolder = [params.datasetFolder 'img/'];

params.outFolder = '../data/';

% minimum distance between two points for which we keep the curve
% fragments.
params.minDist = 3; % in pixels (Eucledean)

params.binSize = 1; % in pixels
% size relative to reference point (adjusted to be a multiple of binSize)
params.relMaxX = 200;
params.relMaxY = 200;
params.relMaxX = floor(params.relMaxX/params.binSize)*params.binSize;
params.relMaxY = floor(params.relMaxY/params.binSize)*params.binSize;
params.relMinX = -params.relMaxX;
params.relMinY = -params.relMaxY;
params.numBins = [(params.relMaxX-params.relMinX)/params.binSize, (params.relMaxY-params.relMinY)/params.binSize];

% get all image names and sizes
files = dir([params.imgsFolder '*jpg']);
params.numImgs = length(files);
params.imgNames = cell(size(files,1),1);
params.imgSizes = zeros(size(files,1),2);
for i=1:size(files,1)
    params.imgNames{i} = files(i).name;
    img = imread([params.imgsFolder files(i).name]);
    params.imgSizes(i,:) = [size(img,1), size(img,2)];
end

% discretization of end point orientation
params.numOrBins = 8;
params.orBinSize = 2*pi/params.numOrBins;


%%
collectCurveFrags(params);

if ~exist('frags')
    load([params.outFolder 'all_frags/frags']);
end
getMeanCurveBetweenInducers(frags, params);

%%
imgNum = 13; % car
% imgNum = 4; % kid
% imgNum = 24; % duck
% imgNum = 1; % plane

if ~exist('frags')
    load([params.outFolder 'all_frags/frags']);
end

demo_completion_single(frags, imgNum, params);