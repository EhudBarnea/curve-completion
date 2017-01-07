function [] = analyzeScaleInvariance(endPointDirection, endPointOr, frags, params)
% compare the center point statistics of curves that belong to inducer
% pairs that are the same up to scale

% endPointDirection - angle between the x axis and the end point
% endPointOr - orientation of the second inducer


visEachScale = false;

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

% keep center point distributions of all scales
densityGridSize = 2^8;
scaleDists = zeros(densityGridSize, densityGridSize, maxScale);
usableScales = false(maxScale, 1);

for s=1:maxScale % loop over scales
%     s
    
    % second inducer
    p2 = [cos(endPointDirection), sin(endPointDirection)] * s;
    or2 = endPointOr;
    x = p2(1);
    y = p2(2);
    
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, false);
    if isUsable
        % rescale returned fragment centers to a cannonical size
        out.fragCenters = (out.fragCenters / s) * 50;
        
        
        % centers point distribution big
        [~,density1,X1,Y1]=kde2d_updated(out.fragCenters, densityGridSize, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.0002);
        % centers point distribution small
        [~,density2,X2,Y2]=kde2d_updated(out.fragCenters, densityGridSize, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.00002);
        if visEachScale
            % show centers point distribution big
            surf(X1,Y1,density1,'EdgeColor','None')
            view(2)
            axis equal
            title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs)]);
            saveas(gcf,[scaleOutFolder 'c_' num2str(s) '_centerD1.png']);
            close all
            
            % show centers point distribution small
            surf(X2,Y2,density2,'EdgeColor','None')
            view(2)
            axis equal
            title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs)]);
            saveas(gcf,[scaleOutFolder 'c_' num2str(s) '_centerD2.png']);
            close all
            
            % show center of completed curve
%             cCenter = c(ceil(size(c,1)/2),:);
%             scatter(cCenter(1),cCenter(2),12,'r','filled')
%             axis equal
%             axis([-200 200 -200 200])
%             title(['num points = ' num2str(size(out.fragCenters,1))]);
%             saveas(gcf,[scaleOutFolder 'c_' num2str(s) '_center.png']);
%             close all
        end
        
        scaleDists(:,:,s) = density2;
        usableScales(s) = true;
    end
end

% calculate distribution difference (from mean) across scales
meanDist = mean(scaleDists(:,:,usableScales),3);
distDiff = abs(scaleDists - repmat(meanDist,1,1,size(scaleDists,3)));
% calculate sum of pixel/bin differences
distSumVar = squeeze(sum(sum(distDiff,2),1)); 
distSumVar(~usableScales) = 0;

% visualize distribution variance across scales
plot(distSumVar)
axis([0 200 0 2])
saveas(gcf,[allScalesOutFolder 'scaleError_' num2str(rad2deg(endPointDirection)) '_' num2str(rad2deg(endPointOr))  '.png']);
close all

end

function [] = createDir(dirName)
% create directory if it doesn't exist

if ~exist(dirName,'dir')
    mkdir(dirName);
end
end

