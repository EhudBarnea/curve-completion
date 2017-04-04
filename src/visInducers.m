function [] = visInducers(p1, or1, p2, or2, onlyOr)
% visualize a pair of inducers, each with its own orientation.

% onlyOr - visualize only the orientation line (without the points)

k = 10; % width of inducer orientation line

firstInducerColor =  [1,0,0];
secondInducerColor = [0,0.7,0];

if nargin < 5
    onlyOr = true;
end

% draw inducer orientations
u = cos(or1);
v = sin(or1);
plot([p1(1),p1(1)-u*k],[p1(2),p1(2)-v*k],'color',firstInducerColor, 'LineWidth', 2);
hold on
u = cos(or2);
v = sin(or2);
% draw inducer points
plot([p2(1),p2(1)-u*k],[p2(2),p2(2)-v*k],'color',secondInducerColor, 'LineWidth', 2);
if ~onlyOr
    hold on
    scatter(p1(1),p1(2),30,firstInducerColor,'filled')
    hold on
    scatter(p2(1),p2(2),30,secondInducerColor,'filled')
end


axis equal
axis([-200 200 -200 200])

end

