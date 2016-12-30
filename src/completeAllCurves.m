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
for x=relMinX:10:relMaxX
    for y=relMinY:10:0 % end at 0, because we mirror curves in y>0
        [x, y]
        
        p2 = [x,y];
        
        % ignore too short curves (there's too many of them)
        if x<=minDist && x>=-minDist && y<=minDist && y>=-minDist
            continue;
        end
        
        for ob=1:params.numOrBins
            
            or2 = getOrFromBin(ob, params.orBinSize);
            [~, isUsable] = completeCurve(p1, or1, p2, or2, frags, params, true);
            if isUsable
                export_fig([params.outFolder '/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '.png']);
            end
            close all;
            
        end
    end
end

end

