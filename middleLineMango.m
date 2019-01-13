function [centroid, headMango, tailMango] = middleLineMango(binMangoImage)

stats = regionprops('table',binMangoImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');

%find centroid of mango with first order image moment
centroid = [stats.Centroid(1,1); stats.Centroid(1,2)];

%find head mango xy coordinate
headMango = [centroid(1) + cosd(stats.Orientation) * 0.5 * stats.MajorAxisLength; centroid(2) - sind(stats.Orientation) * 0.5 * stats.MajorAxisLength];
tailMango = [centroid(1) - cosd(stats.Orientation) * 0.5 * stats.MajorAxisLength; centroid(2) + sind(stats.Orientation) * 0.5 * stats.MajorAxisLength];

