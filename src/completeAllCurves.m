function [] = completeAllCurves(frags, params)
% For each possible inducer pair show all curves and mean curve

minDist = 3;


relMinX = params.relMinX;
relMaxX = params.relMaxX;
relMinY = params.relMinY;
relMaxY = params.relMaxY;


% first inducer is always at point (0,0) looking right
p1 = [0,0];
or1 = 0;

for x=100:10:relMaxX%relMinX:10:relMaxX
    for y=relMinY:10:0 % end at 0, because we mirror curves in y>0
        p2 = [x,y]
        
        % ignore too short curves (there's too many of them)
        if x<=minDist && x>=-minDist && y<=minDist && y>=-minDist
            continue;
        end
        
        for ob=1:params.numOrBins
            
            or2 = getOrFromBin(ob, params.orBinSize);
            [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, true);
            
            % visualize curves and statistics
            if isUsable
                % save figure
                export_fig([params.outFolder '/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '_curves.png']);
                close all;
                
                % show centers point distribution big
                [~,density,X,Y]=kde2d_updated(out.fragCenters, 2^8, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.0002);
                surf(X,Y,density,'EdgeColor','None')
                view(2)
                axis equal
                title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs)]);
                export_fig([params.outFolder '/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '_centerD1.png']);
                close all
                
                % show centers point distribution small
                [~,density,X,Y]=kde2d_updated(out.fragCenters, 2^8, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.00002);
                surf(X,Y,density,'EdgeColor','None')
                view(2)
                axis equal
                title(['num points = ' num2str(size(out.fragCenters,1)) '   num Diff Imgs=' num2str(out.numDiffImgs)]);
                export_fig([params.outFolder '/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '_centerD2.png']);
                close all
                
                % show center of completed curve
                cCenter = c(ceil(size(c,1)/2),:);
                scatter(cCenter(1),cCenter(2),12,'r','filled')
                axis equal
                axis([-200 200 -200 200])
                title(['num points = ' num2str(size(out.fragCenters,1))]);
                export_fig([params.outFolder '/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '_center.png']);
                close all
            end
            close all;
        end
    end
end

end

