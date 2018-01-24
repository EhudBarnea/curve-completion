function [] = demoCompletionSingle(frags, img, params, inducers)
% input two inducers from the user and complete a single curve between them

inputInducers = false;
if nargin < 4
    inputInducers = true;
end

lineWidth = 6;

% load and show image
% close all
imgRGB = img;
hold on
imshow(imgRGB);
set(gcf,'position',get(0,'screensize'));


if inputInducers
    [x,y] = ginput(4);
    
    p1a = [x(1), y(1)];
    p1  = [x(2), y(2)]; % inducer1
    p2a = [x(3), y(3)];
    p2  = [x(4), y(4)]; % inducer2
    
    close all

    % get inducer orientation
    or1 = getOrTwoPts(p1a, p1);
    or2 = getOrTwoPts(p2a, p2);
else
    p1 = [inducers(1), inducers(2)];
    or1 = inducers(3);
    p2 = [inducers(4), inducers(5)];
    or2 = inducers(6);
end



% complete curve
[p1, or1, p2, or2]
[c, isUsable] = completeCurve(p1, or1, p2, or2, frags, params, false);


% draw completion
if isUsable
    figure
    imshow(imgRGB);
    %set(gcf,'position',get(0,'screensize'));

    % draw completion
    % scatter(c(:,1),c(:,2),7,'g','filled');
    hold on
    line(c(:,1),c(:,2),'Color','g', 'LineWidth',lineWidth)
    
    % draw inducers
    hold on
    visInducers(p1, or1, p2, or2, false);
%     line([p1a(1) p1(1)], [p1a(2) p1(2)], 'LineWidth',lineWidth)
%     hold on
%     line([p2a(1) p2(1)], [p2a(2) p2(2)], 'LineWidth',lineWidth)
%     hold on
else
    disp('Completion not usable');
end

end

