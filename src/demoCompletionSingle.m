function [] = demoCompletionSingle(frags, img, params)
% input two inducers from the user and complete a single curve between them


inputInducers = true;

lineWidth = 3;

% load and show image
close all
imgRGB = img;
imshow(imgRGB);
set(gcf,'position',get(0,'screensize'));

% input inducers
if inputInducers
    [x,y] = ginput(4);
    
    p1a = [x(1), y(1)];
    p1  = [x(2), y(2)]; % inducer1
    p2a = [x(3), y(3)];
    p2  = [x(4), y(4)]; % inducer2
else
    p1a = [10, 10];
    p1  = [10, 20]; % inducer1
    p2a = [40, 30];
    p2  = [30, 30]; % inducer2
end

close all

% get inducer orientation
or1 = getOrTwoPts(p1a, p1);
or2 = getOrTwoPts(p2a, p2);

% complete curve
[c, isUsable] = completeCurve(p1, or1, p2, or2, frags, params, false);

% draw completion
if isUsable
    figure
    imshow(imgRGB);
    %set(gcf,'position',get(0,'screensize'));

    % draw completion
    % scatter(c(:,1),c(:,2),7,'g','filled');
    % hold on
    line(c(:,1),c(:,2),'Color','g', 'LineWidth',lineWidth)
    
    % draw inducers
    line([p1a(1) p1(1)], [p1a(2) p1(2)], 'LineWidth',lineWidth)
    hold on
    line([p2a(1) p2(1)], [p2a(2) p2(2)], 'LineWidth',lineWidth)
    hold on
else
    display('Completion not usable');
end

end

