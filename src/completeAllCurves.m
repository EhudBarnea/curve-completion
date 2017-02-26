function [] = completeAllCurves(frags, params)
% For each possible inducer pair show all curves and mean curve

relMinX = params.relMinX;
relMaxX = params.relMaxX;
relMinY = params.relMinY;
relMaxY = params.relMaxY;


% first inducer is always at point (0,0) looking right
p1 = [0,0];
or1 = 0;

xRange = 0:10:relMaxX; %relMinX:10:relMaxX;
parfor xi=1:numel(xRange)
    x=xRange(xi);
    for y=0:-10:relMinY % we mirror curves in y>0
        p2 = [x,y];
        fprintf('%d,%d\n',x,y);
        
        for ob=1:params.numOrBins
            or2 = getOrFromBin(ob, params.orBinSize);
            [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, true);
            
            % visualize curves and statistics
            if isUsable
                % save figure
                saveas(gcf,[params.outFolder 'completions/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '_curves.png']);
                close all;
            end
            close all;
        end
    end
end

end

