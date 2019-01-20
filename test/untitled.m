clc;  % Clear command window.
clear;  % Delete all variables.
close all;  % Close all figure windows except those created by imtool.
workspace;  % Make sure the workspace panel is showing.
fontSize = 20;
x = [0, 560, 0, -60];  
y = [310, 0, -185, 0];
% Append first point to last to close the curve
x = [x, x(1)];
y = [y, y(1)];
plot(x, y, 'r*');
grid on;
set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
set(gcf,'name','Spline Image Analysis Demo','numbertitle','off') 
% Calculate the area within the blue spline curve.
% You do not need to connect the last point back to the first point.
knots = [x; y];
areaOfPolygon = polyarea(x,y);
numberOfPoints = length(x);
% Interpolate with a spline curve and finer spacing.
originalSpacing = 1 : numberOfPoints;
% Make 9 points in between our original points that the user clicked on.
finerSpacing = 1 : 0.1 : numberOfPoints;
% Do the spline interpolation.
splineXY = spline(originalSpacing, knots, finerSpacing);
% Plot the interpolated curve.
hold off;
plot(knots(1, :), knots(2, :), 'ro', 'LineWidth', 2, 'MarkerSize', 16);
hold on;
plot(splineXY(1, :), splineXY(2, :), 'b+-', 'LineWidth', 2, 'MarkerSize', 16);
title('Blue Spline Between Red Knots', 'FontSize', fontSize);
legend('Knots', 'Spline');
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
grid on;
hold off;
% Calculate the area within the blue spline curve.
% You do not need to connect the last point back to the first point.
x = splineXY(1, :);
y = splineXY(2, :);
areaInsideSplineCurve = polyarea(x,y);
% Give the area calculations.
message = sprintf('The area inside the polygon you drew is %.2f.\nThe area inside the blue spline curve is %.2f', ...
  areaOfPolygon, areaInsideSplineCurve);
fprintf(1, '%s', message); % Print to command window.
msgbox(message); % Show user via a popup message box.