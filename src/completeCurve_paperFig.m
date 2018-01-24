function [c, isUsable, out] = completeCurve_paperFig(p1, or1, p2, or2, frags, params, vis)
% complete curve between points p1,p2 with orientations or1,or2

outFolder = '/Users/ehud/Downloads/tmp4/';

% input:
% frags - frags dataset. use fragsSI if params.matchSI=true
% output:
% c - the completed curve (points)
% isUsable - whether the data can be trusted but a good completion
% out - output struct containing.
% out.numFrags - total number of fragments observed between the two inducers
% out.numDiffImgs - number of images from which the used fragments were taken from
% out.fragCenters - center points of all fragments observed between the two inducers

% maxFragsToUse = 300;
maxFragsToUse = inf;
% maxFragsToShow = 10;
maxFragsToShow = 4;
numCurveRepPts = 5;
% numCurveRepPts = 30;


c = [];
out = [];
isUsable = false;

calcCenters = false;
if nargout > 3
    calcCenters = true;
end
    

% vis - visualize completion process
if vis
%     figure
end

% get p2 and or2 in relative to p1 and or1
[endPoint, endPointOr] = transPoints(p2, or2, p1, or1);
% mirror p2
mirrored = false;
if endPoint(2) > 0
    mirrored = true;
    endPoint(2) = -endPoint(2);
    endPointOr = 2*pi - endPointOr;
end

% get relevant fragments
if ~params.matchSI
    endPointFrags = getNearFrags(endPoint, endPointOr, 'rad', frags, params);
else
    endPointFrags = getNearFrags(endPoint, endPointOr, 'si', frags, params);
end

% remove frags from small scales
if params.matchSI
    fragScales = normRows(endPointFrags(:,5:6));
    endPointFrags = endPointFrags(fragScales > params.siMinScale,:);
end

numFrags = size(endPointFrags,1);
numFragsToUse = min(numFrags, maxFragsToUse);
if numFrags < 1
    isUsable = false;
    return;
end

% shuffle curves
% idx = randperm(numFrags);
% endPointFrags = endPointFrags(idx,:);

allRepPts = zeros(numFragsToUse,numCurveRepPts*2); % all curves' representative points
% times 2 because each point is x,y

fragCenters = zeros(numFragsToUse, 2); % center points of all fragments


limitsX = [-30, 110];
limitsY = [-80, 20];
fragsToShow = cell(maxFragsToShow,2);

fragImgs = false(params.numImgs,1); % images with such curves
forFig = [34,18,46,68];
% forFig = [randi(numFragsToUse),randi(numFragsToUse),randi(numFragsToUse),randi(numFragsToUse)];
counter=1;
for i=1:numFragsToUse
    imgNum = endPointFrags(i,1);
    cNum = endPointFrags(i,2); % curve num
    fragP1 = endPointFrags(i,3);
    fragP2 = endPointFrags(i,4);
    annotatorNum = endPointFrags(i,8);
    
    % load curves
    imgName = params.imgNames{imgNum};
    baseName = imgName(1:end-4);
    data = load([params.curvesFolder baseName],'groundTruth');
    curves = data.groundTruth{annotatorNum}; % use one of the annotators
    c = fixCurve(curves{cNum}, params.imgSizes(imgNum,:));
    
    
%     if true||i<=maxFragsToShow
      if any(ismember(forFig,i))
        img = imread([params.imgsFolder imgName]);
        
        fragPts = c;
        p1 = fragP1;
        p2 = fragP2;
        if p1 > p2
%             i
%             annotatorNum
            fragPts = c(end:-1:1,:); % fragment points (flipped)
            p1 = size(fragPts,1)-p1+1;
            p2 = size(fragPts,1)-p2+1;
        end
        
        % flip y axis (so y axis points down)
%         c(:,2) = size(img,2) - c(:,2) + 1;
        dd = fragPts;
        dd(:,2) = size(img,2) - dd(:,2) + 1;
        
        figure
        i
        imshow(img)
        line(dd(:,1),dd(:,2),'Color','b','LineWidth',2);
        hold on
%         visInducers_paperFig(dd(p1,:), or1, dd(p2,:), or2, false);
        visInducersPoint_paperFig(dd(p1,:), or1, dd(p2,:), or2, false);
        hold on
%         line(dd(p1:p2,1),dd(p1:p2,2),'Color','b','LineWidth',2);
        saveas(gcf,[outFolder 'img' num2str(i) '.png']);
%         return
        close all
    end
    
    % transform to canonincal pose
    [fragPts] = getCanonCurve(c, fragP1, fragP2);
    
    
    
    % stretch such that the end point matches endPoint
    if params.relMatchDist || params.matchSI
        fragPts = stretchCurve(fragPts, endPoint);
    end
    
    fragImgs(imgNum) = true;
    
    % get representative points
    repPts = getCurveEquiPoints(fragPts, numCurveRepPts);
    allRepPts(i,:) = reshape(repPts,1,numCurveRepPts*2);
    
    % get frag center
    fragCenters(i,:) = getCurveEquiPoints(fragPts, 1);
    
    % display
%     if vis && i<=maxFragsToShow
    if any(ismember(forFig,i))
        fragsToShow{counter,1} = fragPts;
        fragsToShow{counter,2} = repPts;
        
%         figure
%         line(fragPts(:,1),fragPts(:,2),'Color','k');
%         hold on
%         scatter(repPts(:,1),repPts(:,2),30,'b','filled'); % draw completion points
%         hold on
        
        % ---- show each fragment in its own image
        %         axis equal
        %         axis([-200 200 -200 200])
        %         pause(0.5)
        %         close all
        % ----
        counter = counter+1;
    end
end

% return;

figure;
grid on
line([limitsX(1),limitsX(2)],[0,0],'Color','black','LineWidth',2,'LineStyle','--');
hold on
line([0,0],[limitsY(1),limitsY(2)],'Color','black','LineWidth',2,'LineStyle','--');
hold on
for i=1:maxFragsToShow
    fragPts = fragsToShow{i,1};
    repPts = fragsToShow{i,2};
    line(fragPts(:,1),fragPts(:,2),'Color','k','LineWidth',2);
    hold on
    scatter(repPts(:,1),repPts(:,2),90,'b','filled'); % draw completion points
    hold on
end
hold on
visInducers_paperFig([0, 0], 0, endPoint, endPointOr, false);
axis equal
% axis([-30 130 -100 20])
box on
set(gca,'FontSize',12,'FontWeight','bold','linewidth',2)
xlim(limitsX)
ylim(limitsY)
removeFigureMargin()
% -- fix margin when printing to pdf
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% --

saveas(fig,[outFolder '1_pts.pdf']);
figure

grid on
line([limitsX(1),limitsX(2)],[0,0],'Color','black','LineWidth',2,'LineStyle','--');
hold on
line([0,0],[limitsY(1),limitsY(2)],'Color','black','LineWidth',2,'LineStyle','--');
hold on
scatter(fragCenters(:,1),fragCenters(:,2),10,'b','filled');
hold on
set(gca,'FontSize',12,'FontWeight','bold','linewidth',2)
visInducers_paperFig([0, 0], 0, endPoint, endPointOr, false);
axis equal
xlim(limitsX)
ylim(limitsY)
box on
removeFigureMargin()
% -- fix margin when printing to pdf
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% --
saveas(fig,[outFolder '2_pts.pdf']);
figure

grid on
box on
line([limitsX(1),limitsX(2)],[0,0],'Color','black','LineWidth',2,'LineStyle','--');
hold on
line([0,0],[limitsY(1),limitsY(2)],'Color','black','LineWidth',2,'LineStyle','--');
hold on

% scatter(fragCenters(:,1),fragCenters(:,2),10,'b','filled');
% hold on
% number of different images with such seen curves
numDiffImgs = sum(fragImgs);
isUsable = numFragsToUse>=20;

% prepare output struct
out.fragCenters = fragCenters;
out.numDiffImgs = numDiffImgs;
out.numFrags = numFrags;

% calculate completion curve
c = genCurveMean(allRepPts, numCurveRepPts, params);
% c = genCurveKDE(allRepPts, numCurveRepPts, params);
c = [0,0; c; endPoint];

% split the curve into two parts if not enough samples
if params.extensible && numFrags<params.numMinFrags
    % Get mean center points and orientations
    splitPointIdx = floor(size(c,1))/2;
    splitPoint = c(splitPointIdx,:);
    d=10;
    splitPointBefore = c(splitPointIdx-d,:);
    splitPointAfter = c(splitPointIdx+d,:);
    splitOr1 = getOrTwoPts(splitPointAfter, splitPointBefore);
    splitOr2 = getOrTwoPts(splitPointBefore, splitPointAfter);
    
    params.extensible = false;
    [c1,~,out1] = completeCurve(p1, or1, splitPoint, splitOr1, frags, params, false);
    [c2,~,out2] = completeCurve(splitPoint, splitOr2, p2, or2, frags, params, false);
    c = [c1; c2(2:end,:)];
    numFrags = min([out1.numFrags, out2.numFrags]);
%     scatter(splitPoint(1),splitPoint(2),150,'r')
%     hold on
end

numFrags


if vis
    % draw mean curve
    
    
%     curveColor = [0,0,0.5];
    maxNum = 800;
    cmap = colormap(winter(maxNum));
%     cmap=flipud(cmap);
    curveColor = cmap(min(numFrags,maxNum),:);
    plot(c(:,1), c(:,2), 'Color', curveColor, 'LineWidth' , 3); % draw completed curve
    hold on
    scatter(c(:,1),c(:,2),60,'b','filled'); % draw completion points
    hold on
    visInducers_paperFig([0, 0], 0, endPoint, endPointOr, false);
    set(gca,'FontSize',12,'FontWeight','bold','linewidth',2)
    removeFigureMargin()
    % -- fix margin when printing to pdf
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    % --
%     scatter(0,0,7,'r','filled')
%     hold on
%     scatter(endPoint(1),endPoint(2),7,'b','filled')
    axis equal
%     axis([-30 130 -100 20])
xlim(limitsX)
ylim(limitsY)

saveas(fig,[outFolder '3_pts.pdf']);
%     title(['Num shown curves = ' num2str(min(maxFragsToShow,numFragsToUse)) '   Num used curves = ' num2str(numFragsToUse) '   num Diff Imgs=' num2str(numDiffImgs)])
end

% mirror back curve
if mirrored
    c(:,2) = -c(:,2);
end

% put curve in the coordinate system of p1 and p2
c = transBackPoints(c, p1, or1);

close all
end

function curvePts = genCurveMean(allRepPts, numCurveRepPts, params)
% Generate curve using a set of curves as their mean curve points, i.e.,
% each point in the curve is a mean of all the other corresponding points.

% allRepPts - the set of representative points of the curves to use. This
% is an array of NxM where N is the number of curves and M is
% numCurveRepPts*2 (for x and y).
% numCurveRepPts - the number of representative points of each curve

meanPts = mean(allRepPts,1);
meanPts = reshape(meanPts, numCurveRepPts, 2);

curvePts = meanPts;
end


function curvePts = genCurveKDE(allRepPts, numCurveRepPts, params)
% Generate curve using a set of curves as the maximum of the curve
% distribution calculated using KDE.

% allRepPts - the set of representative points of the curves to use. This
% is an array of NxM where N is the number of curves and M is
% numCurveRepPts*2 (for x and y).
% numCurveRepPts - the number of representative points of each curve


densityGridSize = 2^8;


i=3;
[~,density,X2,Y2]=kde2d_updated(allRepPts(:,i*2:i*2+1), densityGridSize, [params.relMinX, params.relMinY], [params.relMaxX, params.relMaxY], 0.00002);


curvePts = [];
end



function res = normRows(mat)
% calculate the norm of each row in mat
res = sqrt(sum(mat.^2,2));
end



function [] = visInducers_paperFig(p1, or1, p2, or2, onlyOr)
% visualize a pair of inducers, each with its own orientation.

% onlyOr - visualize only the orientation line (without the points)

k = 15; % width of inducer orientation line

firstInducerColor =  [1,0,0];
secondInducerColor = [0,0.7,0];

if nargin < 5
    onlyOr = true;
end

% draw inducer orientations
u = cos(or1);
v = sin(or1);
plot([p1(1),p1(1)-u*k],[p1(2),p1(2)-v*k],'color',firstInducerColor, 'LineWidth', 4);
hold on
u = cos(or2);
v = sin(or2);
% draw inducer points
plot([p2(1),p2(1)-u*k],[p2(2),p2(2)-v*k],'color',secondInducerColor, 'LineWidth', 4);
if ~onlyOr
    hold on
    scatter(p1(1),p1(2),90,firstInducerColor,'filled')
    hold on
    scatter(p2(1),p2(2),90,secondInducerColor,'filled')
end


% axis equal
% axis([-200 200 -200 200])

end

function [] = visInducersPoint_paperFig(p1, or1, p2, or2, onlyOr)
% visualize a pair of inducers, each with its own orientation.

% onlyOr - visualize only the orientation line (without the points)

k = 15; % width of inducer orientation line

firstInducerColor =  [1,0,0];
secondInducerColor = [0,0.7,0];

if nargin < 5
    onlyOr = true;
end

% draw inducer orientations
u = cos(or1);
v = sin(or1);
% plot([p1(1),p1(1)-u*k],[p1(2),p1(2)-v*k],'color',firstInducerColor, 'LineWidth', 2);
% hold on
u = cos(or2);
v = sin(or2);
% draw inducer points
% plot([p2(1),p2(1)-u*k],[p2(2),p2(2)-v*k],'color',secondInducerColor, 'LineWidth', 2);
if ~onlyOr
%     hold on
    scatter(p1(1),p1(2),40,firstInducerColor,'filled')
    hold on
    scatter(p2(1),p2(2),40,secondInducerColor,'filled')
end


% axis equal
% axis([-200 200 -200 200])

end
