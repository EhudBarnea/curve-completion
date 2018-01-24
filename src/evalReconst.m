function [aarc, errors, sortID] = evalReconst(completions, gtCurves, color)
% evaluate reconstruction results over test set. this assumes that the
% first and last points of each curve match those of the ground-truth
% curves

hasColor = false;
if nargin >= 3
    hasColor = true;
end

% ----
% for tests - test only curves of a certain length/disance between inducers
% lenLimits = [100, 120];
% keep = false(length(gtCurves), 1);
% for i=1:length(gtCurves)
%     % get curve
%     gt = gtCurves{i}.pts(gtCurves{i}.p1:gtCurves{i}.p2,:);
%     % Eucledean distance between end points
%     d = sqrt(sum((gt(end,:) - gt(1,:)) .^ 2)); 
%     
%     if d > lenLimits(1) && d < lenLimits(2)
%         keep(i) = true;
%     end
% end
% gtCurves = gtCurves(keep);
% completions = completions(keep);
% ----

% ----
% for tests - test only difficult configurations
% keep = false(length(gtCurves), 1);
% for i=1:length(gtCurves)
%     % get inducers
%     p1 = gtCurves{i}.pts(gtCurves{i}.p1,:);
%     p2 = gtCurves{i}.pts(gtCurves{i}.p2,:);
%     or1 = gtCurves{i}.p1Or;
%     or2 = gtCurves{i}.p2Or;
%     if isDiff(p1, or1, p2, or2)
%         keep(i) = true;
%     end
% end
% gtCurves = gtCurves(keep);
% completions = completions(keep);
% sum(keep)
% ----



% calculate reconstruction errors
errors = zeros(length(gtCurves), 1);
for i=1:length(gtCurves)
    gt = gtCurves{i}.pts(gtCurves{i}.p1:gtCurves{i}.p2,:);
    c = completions{i};

    % get maximal distance between closest points
    [d, dRel] = curve_dist(gt, c);
    errors(i) = dRel;
end

[errors, sortID] = sort(errors); % sort errors
corrects = double(1:length(errors)) / length(errors);

% average of y axis
aarc = 0;
x_points = 0:0.01:1;
for t=x_points
    y=max(corrects(errors<=t));
    if isempty(y)
        y=0;
    end
    aarc = aarc + y/numel(x_points);
end

% error plot
if ~hasColor
    plot(errors, corrects, 'LineWidth',2)
else
    plot(errors, corrects, color, 'LineWidth',2)
end
xlabel('Relative reconstruction error')
ylabel('Ratio of accurately reconstruced curves')
title('Correct reconstructions per permitted error')
xlim([0,1])
ylim([0,1])
grid on
xticks(0:0.2:1)
yticks(0:0.2:1)
set(gca,'FontSize',16)

end