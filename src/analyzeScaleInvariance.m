function [intAngle, meanSTD, meanDiffMu, meanDiffSTD] = analyzeScaleInvariance(endPointDirection, endPointOr, frags, params, visEachScale)
% compare the center point statistics of curves that belong to inducer
% pairs that are the same up to scale

% endPointDirection - angle between the x axis and the end point.
% endPointOr - orientation of the second inducer.
% visEachScale - visualize the center point distribtion at each scale
% (defualt=false).

visScaleSummary = true;

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

% Define maximal scale to examine. At some point we simply don't have
% enough samples.
maxScale = 200;

% Define minimal scale to examine. Very low scales have pretty erronous
% annotations.
minScale = 20;

cannonScale = 50;

removeOutliers = false;

% keep information across all scales
muPts = zeros(maxScale, 2);
scaleSTD = zeros(maxScale, 1);
usableScales = false(maxScale, 1);
fragsPerScale = zeros(maxScale, 1);

% keep all curve centers across all scales
scaleCenters = zeros(0,2);
scaleCentersO = zeros(0,1); % origin scale

for s=31:32%minScale:maxScale % loop over scales
%     s
    
    % second inducer
    p2 = [cos(endPointDirection), sin(endPointDirection)] * s;
    or2 = endPointOr;
    x = p2(1);
    y = p2(2);
    
    [c, isUsable, out, fragCenters] = completeCurve(p1, or1, p2, or2, frags, params, false);
    out.fragCenters = fragCenters;
    if isUsable
        % keep number of fragments
        fragsPerScale(s) = out.numFrags;
        
        % rescale points to a cannonical size
        p2Cannon = ([x,y] / s) * cannonScale;
        out.fragCenters = stretchCurve([out.fragCenters; p2], p2Cannon);
        out.fragCenters = out.fragCenters(1:end-1,:);
        
        % remove outlier centers
        if removeOutliers
            outlierPercent = 0.10;
            mDist = mahal(out.fragCenters,out.fragCenters);
            [~,diffIdx]=sort(-mDist);
            outliersIdx = diffIdx(1:floor(length(diffIdx)*outlierPercent));
            outliers = out.fragCenters(outliersIdx,:);
            out.fragCenters(outliersIdx,:) = [];
        end
        
        % keep centers
        scaleCenters = [scaleCenters; out.fragCenters];
        scaleCentersO = [scaleCentersO; s*ones(size(out.fragCenters,1),1)];
        
        % get mean and STD
        meanCenter = mean(out.fragCenters,1);
        centersSTD = mean(normRows(out.fragCenters-repmat(meanCenter,size(out.fragCenters,1),1)));
        
        if visEachScale
            % visualize center points, inducers, and STD circles
            plot(out.fragCenters(:,1),out.fragCenters(:,2),'b.');
            hold on;
            if removeOutliers
                plot(outliers(:,1),outliers(:,2),'r.');
                hold on;
            end
            scatter(meanCenter(1),meanCenter(2),30,'r','filled');
            hold on;
            viscircles(meanCenter,centersSTD(1), 'LineWidth',1);
            hold on;
            visInducers([0, 0], 0, p2Cannon, or2, false);
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

% p = anova1(scaleCenters',scaleCentersO','off')
% [h,p] = ttest2(scaleCenters(scaleCentersO==31,:),scaleCenters(scaleCentersO==32,:),'Vartype','unequal')
% x=scaleCenters(scaleCentersO==31,2);
% y=scaleCenters(scaleCentersO==32,2);
% [p,h] = ranksum(x,y);


p2 = [cos(endPointDirection), sin(endPointDirection)] * cannonScale;
or2 = endPointOr;

% calculate difference from average mu across scales
meanMu = mean(muPts(usableScales,:),1);
diffMu = normRows(muPts - repmat(meanMu, size(muPts,1), 1));
diffMu(~usableScales) = 0;
meanDiffMu = mean(diffMu(usableScales,:),1);

% calculate std difference from mean std across scales
meanSTD = mean(scaleSTD(usableScales));
diffSTD = abs(scaleSTD - meanSTD);
meanDiffSTD = mean(diffSTD(usableScales));

% calculate line intersection angle representation
intAngle = getLineIntAngle(p2, or2);

if visScaleSummary
    % visualize difference from average mu across scales
    plot(diffMu)
    axis([0 200 0 20])
    saveas(gcf,[allScalesOutFolder num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '_mu.png']);
    close all
    
    % visualize STD at each scales
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
    visInducers([0, 0], 0, p2, or2, false);
    axis([-2 2 -2 2]*cannonScale)
    titleStr = sprintf('%.0f: meanSTD=%.1f, meanDiffMu=%.1f, meanDiffSTD=%.1f',rad2deg(intAngle),meanSTD,meanDiffMu,meanDiffSTD);
    title(titleStr)
    saveas(gcf,[allScalesOutFolder num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '_ind.png']);
    close all
end

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
