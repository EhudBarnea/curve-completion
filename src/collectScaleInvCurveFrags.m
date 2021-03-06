function [] = collectScaleInvCurveFrags(frags, imgNums, params)
% Collect curve fragments in a scale invariant way
% img_nums - list of image numbers from which to draw curves (training
% images)


fragsSIcounts = [];
% Count number of rows for preallocation
[~, fragsSIcounts] = collectAndCount(frags, imgNums, fragsSIcounts, params);
% Collect data in preallocated cell array
[fragsSI, ~] = collectAndCount(frags, imgNums, fragsSIcounts, params);


disp('saving')
if isempty(imgNums)
    save([params.outFolder 'all_frags/fragsSI.mat'],'fragsSI','-v7.3');
else
    save([params.outFolder 'all_frags/fragsSI_train.mat'],'fragsSI','-v7.3');
end
disp('done')
end

function [fragsSI, counts] = collectAndCount(frags, imgNums, fragsSIcounts, params)

onlyCount = false;
if isempty(fragsSIcounts)
    onlyCount = true;
    fragsSI = [];
end

counts = zeros(params.siNumAngBins, params.siNumOrBins);
if ~onlyCount
    fragsSI = cell(params.siNumAngBins, params.siNumOrBins);
    for i = 1:numel(fragsSI)
        fragsSI{i} = zeros(fragsSIcounts(i), 9);
    end
end

% Loop over all the frags array. for fragsSI create a bin for every angular
% distance. put every frag in fragsSI and also keep its angle. when doing
% loopup convert the relative distance to angles
for i = 1:numel(frags)
    if mod(i,1000) == 0
        fprintf('%d/%d, onlyCount=%d\n',i,numel(frags),onlyCount);
    end
    for j = 1:size(frags{i},1)
        c = frags{i}(j,:);
        
        imgNum = frags{i}(j,1);
        cNum = frags{i}(j,2); % curve num
        fragP1 = frags{i}(j,3);
        fragP2 = frags{i}(j,4);
        endPoint = frags{i}(j,5:6);
        endPointOr = frags{i}(j,7);
        annotatorNum = frags{i}(j,8);
        
        % ignore curves from unwanted images
        if ~isempty(imgNums) && ~ismember(imgNum, imgNums)
            continue;
        end
        
        % Get angular location bin of end point.
        % Discretize the [pi,2*pi] region into bins
        angleFromX = 2*pi+atan2(endPoint(2),endPoint(1));
        if angleFromX > 2*pi % correct for endPoint=[x<0,0]
            angleFromX = pi;
        end
        angLocBin = floor((angleFromX-pi)/params.siAngBinSize) + 1;
        if angLocBin > params.siNumAngBins
            angLocBin = params.siNumAngBins;
        end
        
        % get bin of end point orientation
        orBin = getOrBin(endPointOr, params.siOrBinSize, params.siNumOrBins);
        
        c(9) = angleFromX;
        counts(angLocBin,orBin) = counts(angLocBin,orBin) + 1;
        if ~onlyCount
            fragsSI{angLocBin,orBin}(counts(angLocBin,orBin),:) = c;
        end
    end
end

end