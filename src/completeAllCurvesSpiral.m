function [] = completeAllCurvesSpiral(params)
% For each possible inducer pair show all compelted curves, where the
% completion is based on the Euler spiral, calculated  by libcornu and
% saved as a text file.

% read completions file
completions = dlmread([params.outFolder 'completions/libcornu_completions.txt']);


numConfigs = 19*8;
numPts = 200;
for i=1:numConfigs
    
    % get completion points
    ptsStart = 1+(i-1)*(numPts+1)+1;
    ptsEnd = ptsStart+numPts-1;
    pts = completions(ptsStart:ptsEnd,1:2);
    
    % get second inducer
    endPoint = completions(ptsStart-1,3:4);
    endPointOr = completions(ptsStart-1,5);
    
    % flip second inducer orientation
    endPointOr = mod(endPointOr + pi,2*pi);
    
    % draw curve
    maxNum = 800;
    cmap = colormap(winter(maxNum));
    cmap=flipud(cmap);
    curveColor = cmap(maxNum,:);
    plot(pts(:,1), pts(:,2), 'Color', curveColor, 'LineWidth' , 1); % draw completed curve
    hold on;
    visInducers([0, 0], 0, endPoint, endPointOr, false);
    
    if mod(i,19) == 0
        axis equal
        xlim([-250 250]);
        ylim([-300 80]);
        removeFigureMargin();
        
        saveas(gcf,[params.outFolder 'completions/tmp/all_' num2str(endPointOr) '.png']);
        close all;
    end
end



end



% The file 'libcornu_completions.txt' was calculated by adding the
% following function to the file 'test_cornu.c' in libcornu (and executing it from the main function)
% void completeAllCurves(){
%   endpoint e1, e2;
%   cornu_spiral cornu;
% 
%   e1.P[0] = 0;
%   e1.P[1] = 0;
%   e1.theta = 0;
% 
%   float dirRange[] = {180.0010,  190.0000,  200.0000,  210.0000,  220.0000,  230.0000,  240.0000,  250.0000,  260.0000,  270.0000, 280.0000,  290.0000,  300.0000,  310.0000,  320.0000,  330.0000,  340.0000,  350.0000,  360.0000};
%   float orRange[] = {0,    45,    90,   135,   180,   225,   270,   315};
%   float xRange[19];
%   float yRange[19];
%   float scale = 200;
%   int i;
%   for(i=0; i<19; i++){
% 	 xRange[i] = cos(dirRange[i]*M_PI/180)*scale;
% 	 yRange[i] = sin(dirRange[i]*M_PI/180)*scale;
%   }
% 
%   int or2i, p2i;
%   for(or2i=0; or2i<8; or2i++){
%     float or2 = orRange[or2i]*M_PI/180;
%     for(p2i=0; p2i<19; p2i++){
% 
%       e2.P[0] = xRange[p2i];
%       e2.P[1] = yRange[p2i];
%       if(e2.P[0]==0 && e2.P[1]==0)
%         continue;
%       printf("999,999,%f,%f,%f\n", e2.P[0], e2.P[1], e2.theta);
% 
%       // ------------------ this part is copied from the original code
%       cornu_spiral_angles(&cornu, e1.P, e1.theta, e2.P, e2.theta);
% 
%       const int N = 200;
%       int i;
%       double *points, min[2], max[2], scale, D;
%       cairo_surface_t *surface;
%       cairo_t *cr;
%         
% 
%       points = (double *)malloc(sizeof(double) * 2 * (N + 2));
%       cornu_spiral_points(&cornu, 1e-2, N, points);
% 
%       points[2*N] = cornu.P1[0];
%       points[2*N+1] = cornu.P1[1];
%       points[2*N+2] = cornu.P2[0];
%       points[2*N+3] = cornu.P2[1];
% 
%       for(i = 0; i < N; i++) {
%         printf("%lg,%lg\n", points[2*i], points[2*i+1]);
%       }
%       // ------------------
% 
%     }
%   }
% 
% //printf("done\n");
% }
