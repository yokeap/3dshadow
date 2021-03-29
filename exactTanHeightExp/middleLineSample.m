function [centroid, headSample, tailSample] = middleLineSample(binSampleImage)

stats = regionprops('table',binSampleImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');

%find centroid of mango with first order image moment
centroid = [stats.Centroid(1,1); stats.Centroid(1,2)];

%find head mango xy coordinate
headSample = [centroid(1) + cosd(stats.Orientation(1)) * 0.5 * ...
             stats.MajorAxisLength(1); centroid(2) - sind(stats.Orientation(1))...
             * 0.5 * stats.MajorAxisLength(1)];
         
tailSample = [centroid(1) - cosd(stats.Orientation(1)) * 0.5 *...
    stats.MajorAxisLength(1); centroid(2) + sind(stats.Orientation(1)) * ...
    0.5 * stats.MajorAxisLength(1)];

