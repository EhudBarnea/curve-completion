function [] = visInducers(p1, or1, p2, or2, onlyOr)
% visualize a pair of inducers, each with its own orientation.

% onlyOr - visualize only the orientation line (without the points)

k = 15; % width of inducer orientation line

if nargin < 5
    onlyOr = true;
end

% draw inducer orientations
u = cos(or1);
v = sin(or1);
plot([p1(1),p1(1)-u*k],[p1(2),p1(2)-v*k],'color','blue', 'LineWidth', 2);
hold on
u = cos(or2);
v = sin(or2);
% draw inducer points
plot([p2(1),p2(1)-u*k],[p2(2),p2(2)-v*k],'color','blue', 'LineWidth', 2);
if ~onlyOr
    hold on
    scatter(p1(1),p1(2),30,'r','filled')
    hold on
    scatter(p2(1),p2(2),30,'r','filled')
end


axis equal
axis([-200 200 -200 200])

end

