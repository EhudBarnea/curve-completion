function [] = analyzeScaleInvariance(endPointDirection, endPointOr, frags, params, visEachScale)
% compare the center point statistics of curves that belong to inducer
% pairs that are the same up to scale

% endPointDirection - angle between the x axis and the end point.
% endPointOr - orientation of the second inducer.
% visEachScale - visualize the center point distribtion at each scale
% (defualt=false).

if nargin < 5
    visEachScale = false;
end

% first inducer is always at point (0,0) looking right
p1 = [0,0];
or1 = 0;

% create folders for visualization
createDir([params.outFolder 'stats/scale_analysis/']);
allScalesOutFolder = [params.outFolder 'stats/scale_analysis/all_scales/'];
createDir(allScalesOutFolder);
if visEachScale
    scaleOutFolder = [params.outFolder 'stats/scale_analysis/scale_' num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr)) '/'];
    createDir(scaleOutFolder);
end

% define maximal scale to examine
maxScale = 200;

cannonScale = 50;

removeOutliers = false;

% keep information across all scales
muPts = zeros(maxScale, 2);
scaleSTD = zeros(maxScale, 1);
usableScales = false(maxScale, 1);
fragsPerScale = zeros(maxScale, 1);

for s=1:maxScale % loop over scales
%     s
    
    % second inducer
    p2 = [cos(endPointDirection), sin(endPointDirection)] * s;
    or2 = endPointOr;
    x = p2(1);
    y = p2(2);
    
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, false);
    if isUsable
        % keep number of fragments
        fragsPerScale(s) = out.numFrags;
        
        % rescale points to a cannonical size
        out.fragCenters = (out.fragCenters / s) * cannonScale;
        p2Cannon = ([x,y] / s) * cannonScale;
        
        % remove outlier centers
        if removeOutliers
            outlierPercent = 0.10;
            mDist = mahal(out.fragCenters,out.fragCenters);
            [~,diffIdx]=sort(-mDist);
            outliersIdx = diffIdx(1:floor(length(diffIdx)*outlierPercent));
            outliers = out.fragCenters(outliersIdx,:);
            out.fragCenters(outliersIdx,:) = [];
        end
        
        % fit gaussian (get std along principal directions)
        meanCenter = mean(out.fragCenters,1);
        centersSTD = mean(normRows(out.fragCenters-repmat(meanCenter,size(out.fragCenters,1),1)));
        
        if visEachScale
            % visualize center points, inducers, and STD circles
            plot(out.fragCenters(:,1),out.fragCenters(:,2),'b.');
            hold on;
            plot(outliers(:,1),outliers(:,2),'r.');
            hold on;
            scatter(meanCenter(1),meanCenter(2),30,'r','filled');
            hold on;
            viscircles(meanCenter,centersSTD(1), 'LineWidth',1);
            hold on;
            visInducers([0, 0], 0, p2Cannon, or2);
            axis([-20,100,-70,50]);
            title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs) '   STD1=' num2str(centersSTD(1))]);
            saveas(gcf,[scaleOutFolder 'c_' sprintf('%03d',s) '_pts.png']);
            close all;
        end
        
        muPts(s,:) = meanCenter;
        scaleSTD(s) = centersSTD(1);
        usableScales(s) = true;
    end
end


% calculate difference from average mu across scales
meanMu = mean(muPts(usableScales,:),1);
diffMu = normRows(muPts - repmat(meanMu, size(muPts,1), 1));
diffMu(~usableScales) = 0;

% visualize difference from average mu across scales
plot(diffMu)
axis([0 200 0 20])
saveas(gcf,[allScalesOutFolder num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '_mu.png']);
close all

% visualize Gaussian STD at each scales
plot(scaleSTD)
axis([0 200 0 30])
saveas(gcf,[allScalesOutFolder num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '_std.png']);
close all

% visualize number of fragments at each scale
plot(fragsPerScale)
saveas(gcf,[allScalesOutFolder num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '_nFrags.png']);
close all

% visualize inducers
scatter(meanMu(1),meanMu(2),30,'r','filled');
hold on
p2 = [cos(endPointDirection), sin(endPointDirection)] * cannonScale;
or2 = endPointOr;
visInducers([0, 0], 0, p2, or2);
axis([-2 2 -2 2]*cannonScale)
saveas(gcf,[allScalesOutFolder num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '_ind.png']);
close all

end

function [] = createDir(dirName)
% create directory if it doesn't exist

if ~exist(dirName,'dir')
    mkdir(dirName);
end
end


function res = normRows(mat)
    % calculate the norm of each row in mat
    res = sqrt(sum(mat.^2,2));
end
