function [splineXY] = splineEstimate(XY)

    XY = XY';
    % Append first point to last to close the curve
    %XY = [XY(1,:), XY(1,1); XY(2,:), XY(2,1)];
    
    knots = [XY(1,:); XY(2,:)];
    numberOfPoints = length(XY(1,:));
    
    % Interpolate with a spline curve and finer spacing.
    originalSpacing = 1 : numberOfPoints;
    % Make 9 points in between our original points that the user clicked on.
    finerSpacing = 1 : 0.1 : numberOfPoints;
    % Do the spline interpolation.
    splineXY = spline(originalSpacing, knots, finerSpacing);
    splineXY = splineXY';
% 
%     % Calculate the area within the spline curve.
%     x = splineXY(1, :);
%     y = splineXY(2, :);
%     areaInsideSplineCurve = polyarea(x,y);
    
end

