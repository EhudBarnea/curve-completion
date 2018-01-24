function [c] = completeCurveCS(p1, or1, p2, or2)
% complete curve between points p1,p2 with orientations or1,or2 with a
% cubic spline

vis = true;

x = [p1(1), p2(1)];
y = [p1(2), p2(2)];

if x(1) ~= x(2)
    or2t = mod(or2 + pi,2*pi);
    cs = spline(x,[or1 y or2t]);
    xx = linspace(p1(1), p2(1),101);
    yy = ppval(cs,xx);
else
    yy = min([y(1), y(2)]):0.1:max([y(1), y(2)]);
    xx = ones(1, numel(yy)) * x(1);
end
c = [xx', yy'];

if vis
    plot(c(:,1), c(:,2), 'Color', 'b', 'LineWidth' , 1); % draw completed curve
    hold on
    visInducers(p1, or1, p2, or2, false);
    axis equal
    axis([-200 200 -200 200])
end

% plot(xx,yy);
% xlim([-2,2])
% ylim([-2,2])


end
