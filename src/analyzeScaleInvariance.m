function [] = analyzeScaleInvariance(endPointDirection, endPointOr, frags, params)
% compare the center point statistics of curves that belong to inducer
% pairs that are the same up to scale

% endPointDirection - direction from the first inducer at (0,0) to the second inducer
% endPointOr - orientation of the second inducer


or2 = endPointOr;

scaleOutFolder = [params.outFolder 'stats/scale_' num2str(endPointDirection(1)) '_' num2str(endPointDirection(2)) '_' num2str(or2) '/'];
mkdir(scaleOutFolder);

% first inducer is always at point (0,0) looking right
p1 = [0,0];
or1 = 0;

for s=1:100 % loop over scales
    s
    
    p2 = round(endPointDirection * s);
    x = p2(1);
    y = p2(2);
    
    [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, false);
    if isUsable
        % rescale returned fragment centers to a cannonical size
        out.fragCenters = (out.fragCenters / s) * 50;
        
        % show centers point distribution big
        [~,density,X,Y]=kde2d_updated(out.fragCenters, 2^8, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.0002);
        surf(X,Y,density,'EdgeColor','None')
        view(2)
        axis equal
        title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs)]);
        export_fig([scaleOutFolder 'c_' num2str(x) '_' num2str(y) '_' num2str(or2) '_centerD1.png']);
        close all
        
        % show centers point distribution small
        [~,density,X,Y]=kde2d_updated(out.fragCenters, 2^8, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.00002);
        surf(X,Y,density,'EdgeColor','None')
        view(2)
        axis equal
        title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs)]);
        export_fig([scaleOutFolder 'c_' num2str(x) '_' num2str(y) '_' num2str(or2) '_centerD2.png']);
        close all
        
        % show center of completed curve
        cCenter = c(ceil(size(c,1)/2),:);
        scatter(cCenter(1),cCenter(2),12,'r','filled')
        axis equal
        axis([-200 200 -200 200])
        title(['num points = ' num2str(size(out.fragCenters,1))]);
        export_fig([scaleOutFolder 'c_' num2str(x) '_' num2str(y) '_' num2str(or2) '_center.png']);
        close all
    end
end

