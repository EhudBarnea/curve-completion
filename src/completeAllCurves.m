function [] = completeAllCurves(frags, params)
% For each possible inducer pair show all mean curve

scaleInvariant = true;

relMinX = params.relMinX;
relMaxX = params.relMaxX;
relMinY = params.relMinY;
relMaxY = params.relMaxY;


% first inducer is always at point (0,0) looking right
p1 = [0,0];
or1 = 0;

if ~scaleInvariant
    skip = 50;
    xRange = relMinX:skip:relMaxX-skip; %relMinX:10:relMaxX;
    yRange = 0:-skip:relMinY+skip; % we mirror curves in y>0
    [m1, m2] = meshgrid(xRange, yRange);
    xyRange = [m1(:), m2(:)];
    
    orRange = deg2rad(0:45:315);
else
%     scale = 80;
    scale = 200;
    % loop over different end point directions and end point orientations. The range variables are used since
    % parfor requires consecutive integers.
    dirRange = deg2rad([180.001, 190:10:360])';
    orRange = deg2rad(0:45:315);
    
    % get points from direction
    xyRange = [cos(dirRange), sin(dirRange)] * scale;
end


for or2i=1:numel(orRange)
% parfor or2i=numel(orRange)
    or2 = orRange(or2i);
    for p2i=1:size(xyRange,1)
        
        
        x=xyRange(p2i,1);
        y=xyRange(p2i,2);
        if x==0 && y==0
            continue;
        end
        p2 = [x,y];
        fprintf('%f,%f\n',x,y);
        
        
        [c, isUsable, out] = completeCurve(p1, or1, p2, or2, frags, params, true);
        
        hold on
        title('');
        axis equal
%         xlim([-220 220]); 
%         ylim([-230 50]); 
        xlim([-250 250]);
        ylim([-300 80]);
        saveas(gcf,[params.outFolder 'completions/tmp/all_' num2str(or2) '.png']);
        % visualize curves and statistics
        if false
            if isUsable
                % save figure
                saveas(gcf,[params.outFolder 'completions/tmp/c_' num2str(x) '_' num2str(y) '_' num2str(ob) '_curves.png']);
                close all;
            end
            close all;
        end
    end
    
    close all;
end

disp('done');
end

