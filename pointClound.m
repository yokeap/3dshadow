clc;  % Clear command window.
clear;  % Delete all variables.
close all;  % Close all figure windows except those created by imtool.
workspace;  % Make sure the workspace panel is showing.
img = imread('DSC_0234.JPG');
[binMangoImage, maskedMangoRGBImage] = segmentImage(img);
[binMangoShadowImage, maskedMangoShaowRGBImage] = segmentShadow(img);


binMangoImage = bwareafilt(binMangoImage, 1);                              %filter by area threshold
binMangoImage = imfill(binMangoImage,'holes');                              %fill holes
binMangoShadowImage = imfill(binMangoShadowImage,'holes');                 %fill holes
%binMangoShadowImage = bwareaopen(binMangoShadowImage, 500);                  %filter by area threshold

% for i = 1:50
%     binaryImage = medfilt2(binaryImage,[50 50]);
% end

%binaryImage = medfilt2(binaryImage,[100 100]);

%apply median filtering with 10x10 kernel for smoothing edge
binMangoImage = medfilt2(binMangoImage,[10 10]);
binMangoShadowImage = medfilt2(binMangoShadowImage, [10 10]);
binMangoShadowImage = bwareaopen(binMangoShadowImage, 500);

%Overlay of binary image of mango and shadow
rgbOverlay = imfuse(binMangoImage, binMangoShadowImage);

%Extracting of edges of mango and shadow
binMangoEdgeImage = edge(binMangoImage,'Canny');   
binMangoShadowEdgeImage = edge(binMangoShadowImage, 'Canny');

[MangoEdgesCol,MangoEdgesRow] = find(binMangoEdgeImage);                   %find the non-zero element (return [x,y])
[MangoShadowEdgesCol,MangoShadowEdgesRow] = find(binMangoShadowEdgeImage);                   %find the non-zero element (return [x,y])

%extract edges with maximum amplitude and lowest amplitude
mangoEdges = extractEdges(binMangoImage);
%extract shadow edges based on positon of mango edges
mangoShadowEdges = extractShadowEdges(binMangoShadowImage, mangoEdges(:,1));

%compute length of shadow
mangoShadowLength = computeShadowLength(mangoShadowEdges);

% plot mango edges
% plot(MangoEdgesRow,MangoEdgesCol,'.');
% hold on; 
% plot(MangoShadowEdgesRow,MangoShadowEdgesCol,'.');


% plot(mangoEdges(:,1), mangoEdges(:,2), '.');
% hold on
% plot(mangoEdges(:,1), mangoEdges(:,3), '.');
% 
% [centroid, headMango, tailMango] = middleLineMango(binMangoImage);
% 
% hold on
% plot([headMango(1), tailMango(1)], [headMango(2), tailMango(2)], '-');

%trim vector lenght
n = min(length(mangoEdges(:,1)), length(mangoShadowEdges(:,1)));
middleHeight = round(min(mangoShadowEdges(1:n,3)) / 2);

% mangoShadowLength(1:n,2) = round(mangoShadowLength(1:n,2) * 10);
% middleHeight = min(mangoShadowLength(1:n,2)) / 2;

middleClound(:,1) = mangoEdges(1:n,1);
middleClound(:,2) = mangoEdges(1:n,2);
middleClound(:,3) = mangoEdges(1:n,3);
middleClound(:,4) = 0;

%estimte centroid line of mango with linear approximation by 2 points
%which was estimated from first order image moments centroid
[centroid, headMango, tailMango] = middleLineMango(binMangoImage);
coefficients = polyfit([tailMango(1), headMango(1)], [tailMango(2), headMango(2)], 1);
polyMiddleLine = polyval(coefficients, mangoEdges(1:n,1));

upperClound(:,1) = mangoEdges(1:n,1);
upperClound(:,2) = polyMiddleLine(1:n,1);
%upperClound(:,3) = round((mangoShadowEdges(1:n,3) / 2) + middleHeight);
upperClound(:,3) = round(mangoShadowLength(1:n,2) / 2);

lowestClound(:,1) = mangoEdges(1:n,1);
lowestClound(:,2) = polyMiddleLine(1:n,1);
lowestClound(:,3) = upperClound(:,3) * -1;
%lowestClound(:,3) = abs(round(middleHeight - (upperClound(:,3) - middleHeight)));
%lowestClound(:,3) = abs((mangoShadowEdges(1:n,3) / 2) - middleHeight);

scatter3(middleClound(:,1), middleClound(:,2), middleClound(:,4), 'filled')
hold on
scatter3(middleClound(:,1), middleClound(:,3), middleClound(:,4), 'filled')
scatter3(upperClound(:,1), upperClound(:,2), upperClound(:,3), 'filled')
scatter3(lowestClound(:,1), lowestClound(:,2), lowestClound(:,3), 'filled')

%stroring all of 3D data in 7x1 array
%%% [ slicing amount(x co.); 
%%%   middle height of object; 
%%%   central line of object in x co.; 
%%%   low side object in y co.;
%%%   high side object in y co.;
%%%   height side of object in z co.;
%%%%  low side of object in z co.;
%%%% ]
ptClound = [mangoEdges(1:n,1),middleClound(:,4), polyMiddleLine(1:n,1), mangoEdges(1:n,2), mangoEdges(1:n,3), upperClound(:,3), lowestClound(:,3)];

% %New figure establishing
% figure
% %100th slicing plot
% plot(ptClound(100,4:5), ptClound(100,2), '*')
% hold on
% plot(ptClound(100,3), ptClound(100,6), '*') 
% plot(ptClound(100,3), ptClound(100,7), '*')
% %center of slice
% plot(ptClound(100,3), ptClound(100,2), '*')

plotCrossSection(ptClound, 200);

%Slicing sampling with counter clockwise
% slicing(:,1) = [ptClound(100,4); ptClound(100,3); ptClound(100,5); ptClound(100,3)];
% slicing(:,2) = [ptClound(100,2); ptClound(100,6); ptClound(100,2); ptClound(100,7)];

figure
hold on
%half slide along in z axis storing.
leftSliceX = [];
leftSliceY = [];
leftSliceZ = [];
rightSliceX = [];
rightSliceY = [];
rightSliceZ = [];
for i = 1:n
    %left half slice
    halfSlice_L(:,1,i) = [ptClound(i,3); ptClound(i,4); ptClound(i,3)];
    halfSlice_L(:,2,i) = [ptClound(i,6); ptClound(i,2); ptClound(i,7)];
    %right half slice
    halfSlice_R(:,1,i) = [ptClound(i,3); ptClound(i,5); ptClound(i,3)];
    halfSlice_R(:,2,i) = [ptClound(i,6); ptClound(i,2); ptClound(i,7)];
    
    %spline interpolation
    sliceModel_L(:,2:3,i) = splineEstimate(halfSlice_L(:,:,i));
    sliceModel_R(:,2:3,i) = splineEstimate(halfSlice_R(:,:,i));
    %append x coordinate to column 1
    sliceModel_L(:,1,i) = mangoEdges(i,1);
    sliceModel_R(:,1,i) = mangoEdges(i,1);
    
    %scatter3(sliceModel_R(:,1,i), sliceModel_R(:,2,i), sliceModel_R(:,3,i), 2)
    leftSliceX = [leftSliceX; sliceModel_R(:,1,i)];
    leftSliceY = [leftSliceY; sliceModel_R(:,2,i)];
    leftSliceZ = [leftSliceZ; sliceModel_R(:,3,i)];
    
    %scatter3(sliceModel_R(:,1,i), sliceModel_R(:,2,i), sliceModel_R(:,3,i), 2)
    rightSliceX = [rightSliceX; sliceModel_L(:,1,i)];
    rightSliceY = [rightSliceY; sliceModel_L(:,2,i)];
    rightSliceZ = [rightSliceZ; sliceModel_L(:,3,i)];
end

slice = 200;

figure
plot(sliceModel_L(:,2,slice), sliceModel_L(:,3,slice), 'r-');
hold on;
plot(sliceModel_R(:,2,slice), sliceModel_R(:,3,slice), 'b-');

figure
plot3(leftSliceX,leftSliceY,leftSliceZ, 'ro');
hold on
plot3(rightSliceX,rightSliceY,rightSliceZ, 'bo');