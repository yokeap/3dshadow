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
middleHeight = max(mangoShadowEdges(1:n,3)) / 2;

middleClound(:,1) = mangoEdges(1:n,1);
middleClound(:,2) = mangoEdges(1:n,2);
middleClound(:,3) = mangoEdges(1:n,3);
middleClound(:,4) = middleHeight;

%estimte centroid line of mango with linear approximation by 2 points
%which was estimated from first order image moments centroid
[centroid, headMango, tailMango] = middleLineMango(binMangoImage);
coefficients = polyfit([tailMango(1), headMango(1)], [tailMango(2), headMango(2)], 1);
polyMiddleLine = polyval(coefficients, mangoEdges(1:n,1));

upperClound(:,1) = mangoEdges(1:n,1);
upperClound(:,2) = polyMiddleLine(1:n,1);
upperClound(:,3) = round((mangoEdges(1:n,3) / 2) + middleHeight);

lowestClound(:,1) = mangoEdges(1:n,1);
lowestClound(:,2) = polyMiddleLine(1:n,1);
lowestClound(:,3) = abs(round((mangoEdges(1:n,3) / 2) - middleHeight));

scatter3(middleClound(:,1), middleClound(:,2), middleClound(:,4), 'filled')
hold on
scatter3(middleClound(:,1), middleClound(:,3), middleClound(:,4), 'filled')
scatter3(upperClound(:,1), upperClound(:,2), upperClound(:,3), 'filled')
scatter3(lowestClound(:,1), lowestClound(:,2), lowestClound(:,3), 'filled')

ptClound = [mangoEdges(1:n,1), middleHeight, polyMiddleLine(1:n,1), mangoEdges(1:n,2), mangoEdges(1:n,3), round((mangoEdges(1:n,3) / 2) + middleHeight), abs(round((mangoEdges(1:n,3) / 2) - middleHeight))];

%100th slicing plot
plot(ptClound(100,4:5), ptClound(100,2), '*')
hold on
plot(ptClound(100,3), ptClound(100,6), '*')
plot(ptClound(100,3), ptClound(100,7), '*')
%center of slice
plot(ptClound(100,3), ptClound(100,2), '*')

%example slicing
slicing(:,1) = [ptClound(100,4); ptClound(100,5); ptClound(100,3); ptClound(100,3)];
slicing(:,2) = [ptClound(100,2); ptClound(100,2); ptClound(100,6); ptClound(100,7)];

