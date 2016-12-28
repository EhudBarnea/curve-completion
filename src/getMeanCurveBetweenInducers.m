function [] = getMeanCurveBetweenInducers(frags, params)
% For each two inducers show all curves and mean curve

numCurveRepPts = 5;
vis = true;
minDist = 3;


relMinX = params.relMinX;
relMaxX = params.relMaxX;
relMinY = params.relMinY;
relMaxY = params.relMaxY;
numOrBins = params.numOrBins;
binSize = params.binSize;
numBins = params.numBins;
numImgs = params.numImgs;
imgNames = params.imgNames;
curvesFolder = params.curvesFolder;
imgSizes = params.imgSizes;
outFolder = params.outFolder;

for x=relMinX:10:relMaxX
    for y=relMinY:10:0 % end at 0, because we mirror curves in y>0
        [x, y]
        
        % ignore too short curves (there's too many of them)
        if x<=minDist && x>=-minDist && y<=minDist && y>=-minDist
            continue;
        end
        
        for ob=1:numOrBins
            endPoint = [x,y];
            endPointBin = floor((endPoint - [relMinX, relMinY])/binSize) + 1;
            endPointBin(endPointBin>numBins) = numBins(1);
            endPointFrags = frags{endPointBin(1), endPointBin(2), ob};
            
            numFrags = size(endPointFrags,1);
            if numFrags < 1
                continue;
            end
            
            % shuffle curves
            idx = randperm(numFrags);
            endPointFrags = endPointFrags(idx,:);
            
            allRepPts = zeros(numFrags,numCurveRepPts*2); % all curves' representative points
            % times 2 because each point is x,y
            
            maxFragsToShow = 10;
            fragImgs = false(numImgs,1); % images with such curves
            for i=1:numFrags
                imgNum = endPointFrags(i,1);
                cNum = endPointFrags(i,2); % curve num
                p1 = endPointFrags(i,3);
                p2 = endPointFrags(i,4);
                
                % load curves
                imgName = imgNames{imgNum};
                baseName = imgName(1:end-4);
                data = load([curvesFolder baseName],'groundTruth');
                curves = data.groundTruth{1}; % use set 1 (each set is a different annotator)
                c = fixCurve(curves{cNum}, imgSizes(imgNum,:));
                
                % transform to canonincal pose
                [fragPts] = getCanonCurve(c, p1, p2);
                
                fragImgs(imgNum) = true;
                
                % get representative points
                repPts = getCurveEquiPoints(fragPts, numCurveRepPts);
                allRepPts(i,:) = reshape(repPts,1,numCurveRepPts*2);
                
                % display
                if vis && i<=maxFragsToShow
                    line(fragPts(:,1),fragPts(:,2));
                    hold on
                end
            end
            
            % number of different images with such seen curves
            numDiffImgs = sum(fragImgs);
            
            % calculate mean curve
            %         coeff = pca(allRepPts);
            meanPts = mean(allRepPts,1);
            meanPts = reshape(meanPts, numCurveRepPts, 2);
            meanPts = [[0,0]; meanPts; endPoint];
            
            meanCurve = meanPts;
            %             save(['curves_between_inducers/' num2str(endPoint(1)) '_' num2str(endPoint(2)) '_' num2str(orBin) '.mat'],'meanCurve');
            
            
            if vis && numDiffImgs>2
                % draw mean curve
                scatter(meanPts(:,1),meanPts(:,2),7,'r','filled');
                line(meanPts(:,1),meanPts(:,2),'Color','r');
                
                hold on
                scatter(0,0,7,'r','filled')
                hold on
                scatter(endPoint(1),endPoint(2),7,'r','filled')
                axis equal
                axis([-200 200 -200 200])
                title(['Num shown curves = ' num2str(min(40,numFrags)) '   num Diff Imgs=' num2str(numDiffImgs)])
                
                export_fig([outFolder '/curves_between_inducers/c_' num2str(endPoint(1)) '_' num2str(endPoint(2)) '_' num2str(ob) '.png']);
            end
            close all;
        end
    end
end

end

