function [resCurves] = completeAllCurves(frags, params)
% For each possible inducer pair show all mean curve


% first inducer is always at point (0,0) looking right
p1 = [0,0];
or1 = 0;


% relMinX = params.relMinX;
% relMaxX = params.relMaxX;
% relMinY = params.relMinY;
% relMaxY = params.relMaxY;
% skip = 50;
% xRange = relMinX:skip:relMaxX-skip; %relMinX:10:relMaxX;
% yRange = 0:-skip:relMinY+skip; % we mirror curves in y>0
% [m1, m2] = meshgrid(xRange, yRange);
% xyRange = [m1(:), m2(:)];
% orRange = deg2rad(0:45:315);



if ~params.matchSI
    scale=80;
else
    scale = 200;
end
% loop over different end point directions and end point orientations. The range variables are used since
% parfor requires consecutive integers.
dirRange = deg2rad([180.001, 190:10:360])';
orRange = deg2rad(0:45:315);

% get points from direction
xyRange = [cos(dirRange), sin(dirRange)] * scale;

% keep all generated curves and their number of frags
resCurves = {};

for or2i=1:numel(orRange)
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
        resCurves{length(resCurves)+1}.pts = c;
        resCurves{length(resCurves)}.numFrags = out.numFrags;
        
        hold on
        title('');
        axis equal
        if ~params.matchSI
            k = 0.4; % for non scale invariant
        else
            k = 1;
        end
        xlim([-250 250]*k);
        ylim([-300 80]*k);

        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        box on
        set(gca,'FontSize',25,'FontWeight','bold','linewidth',2)
        fig = gcf;
        if p2i==1%p2i==size(xyRange,1)
%             removeFigureMargin();
%             set(gca,'LooseInset',get(gca,'TightInset'))
            % -- fix margin when printing to pdf
            fig.PaperPositionMode = 'auto';
            fig_pos = fig.PaperPosition;
            fig.PaperSize = [fig_pos(3) fig_pos(4)];
            % --
            
        end
        model_out_folder = 'tmp';
%         model_out_folder = 'mean_curve_100';
%         model_out_folder = 'mean_curve_si';
%         model_out_folder = 'mean_curve_si_ext';
%         model_out_folder = 'mean_curve_si_ext_with_point';
        outFilename = [params.outFolder 'completions/' model_out_folder '/' num2str(or2i)];
        %         saveas(gcf,[outFilename '.png']);
        %         saveas(gcf,[outFilename '.eps'], 'epsc2');
        print(fig, '-dpdf', '-r600', [outFilename '.pdf'])
%         export_fig([outFilename '.pdf']);
        
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

