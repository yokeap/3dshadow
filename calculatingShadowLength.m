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
%Shadow height calculation
mangoEdges = extractEdges(binMangoImage);
mangoShadowEdges = extractEdges(binMangoShadowImage);

% plot mango edges
% plot(MangoEdgesRow,MangoEdgesCol,'.');
% hold on; 
% plot(MangoShadowEdgesRow,MangoShadowEdgesCol,'.');


plot(mangoEdges(:,1), mangoEdges(:,2), '.');
hold on
plot(mangoEdges(:,1), mangoEdges(:,3), '.');

mangoShadowEdgesLength = 
