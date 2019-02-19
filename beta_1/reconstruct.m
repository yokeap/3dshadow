close all;
clear;
clc;

homographyMatrix = ...
[ -0.000001171027520  -0.000569119373972   0.402795558322037
  -0.000571860777008  -0.000007045512653   0.915276425974051
  -0.000000000072476  -0.000000017684098  -0.004914395139275];

% virtual light source position
virLightPos = [-2.298030767150990e+02,1.469862104150302e+02,4.943304818724096e+02];

% Reference tower position in image coordinate.
PA = [1.591724358974359e+03,1.427596153846153e+03;
      6.980448717948707e+02,9.553269230769220e+02;
      1,1]';
  
PB = [3.876211538461539e+03,3.995980769230770e+03;
      6.980448717948707e+02,9.553269230769220e+02;
      1,1]';

PC = [3.876211538461539e+03,3.995980769230770e+03;
      3.177711538461538e+03,3.705583333333333e+03;
      1,1]';

PD = [1.565108974358974e+03,1.409852564102564e+03;
      3.177711538461538e+03,3.705583333333333e+03;
      1,1]';

% Refernce tower position in World coordinate.
HA = [-0.344940596426928,-0.030344488603148,0;   % bottom of blue tower corner
     28.500421768196030,-18.857746363838732,0;  % shadow ray casted from top of blue tower corner
     -0.344940596426928,-0.030344488603148,52]; % height of top of blue tower corner 
  
HB = [0.197484825228561,2.648525869040325e+02,0;
     29.111778283175088,2.802209834378628e+02,0;
     0.197484825228561,2.648525869040325e+02,52];
 
HC = [2.833449765829964e+02,2.660276503092134e+02,0;
     3.429544363039715e+02,2.799133485893511e+02,0;
     2.833449765829964e+02,2.660276503092134e+02,52];
 
HD = [2.828118364137633e+02,1.035175223803500,0;
     3.437996501867177e+02,-16.227331116713668,0;
     2.828118364137633e+02,1.035175223803500,52];
 
% azimuth triangle (XY) or top view in radian
azimuth_a = 1.3456;
azimuth_b = 1.2907;
azimuth_c = 0.5053;     % this angle is used for ray casting in top view
 
% elevation triangle (XZ) or side view in radian
elevation_a = 0.7113;
elevation_b = 2.0523;
elevation_c = 0.3780;

 % load intrinsic camera calibrated parameters
load('../calibration/cameraParams.mat');

% load distorted image
imOrig = imread('../SampleImages/background.JPG');
% get undistorting image and offset origin
[im, newOrigin] = undistortImage(imOrig, cameraParams, 'OutputView', 'full');
figure; imshow(im); title('Undistorted Image');
axis image;
axis ij ;

% load background image.
imgBg = imread('../SampleImages/background.JPG');
% load sample image.
imgSample = imread('../SampleImages/DSC_0326.JPG');

% seperate object by subtracts sample image and background.
diffImage = imsubtract(imgBg, imgSample);
% convert to binary and filtering image.
[diffImageBW, maskedRGBImage] = createMaskBW(diffImage);
diffImageBW = bwareafilt(diffImageBW, 1);
diffImageBW = imfill(diffImageBW,'holes');
imshow(diffImageBW);

% create blank blackgroumd image (background color is according to color of object).
grennScreenImage=zeros(size(imgSample,1), size(imgSample,2), 3,'uint8');
grennScreenImage(:,:,2) = 255;

% masking the input image from masked image that come from subtracted image.
% *note: matlab direct accessing pixel is (y,x), ie img(y,x).
for i=1:size(imgSample,2)
    for j=1:size(imgSample,1)
        if(diffImageBW(j,i) > 0)
            grennScreenImage(j,i,:) = imgSample(j,i,:);
        end      
    end
end
 
figure;
imshow(grennScreenImage);

% Sample segmentation,
%---------------------------------------------------------------------
% color thresholding (HSV model)
[binSegmentImage, MaskedSegmentRGBImage] = appleMask(grennScreenImage);        
% filtering for only 1 largest object
binSegmentImage = bwareafilt(binSegmentImage, 1);             
% fill hole
binSegmentImage = imfill(binSegmentImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentImage  = medfilt2(binSegmentImage ,[3 3]);
figure; imshow(binSegmentImage);

% 
% %Extracting of edges of sample and shadow
% binSampleEdgeImage = edge(binSegmentImage,'Canny');   
% binSegmentShadowEdgeImage = edge(binSegmentShadowImage, 'Canny');
% figure; imshow(binSampleEdgeImage);
% figure; imshow(binSegmentShadowEdgeImage);
% 
% %find the non-zero element (return [x,y])
% [SampleEdgesCol, SampleEdgesRow] = find(binSampleEdgeImage);
% [SampleShadowEdgesCol, SampleShadowEdgesRow] = find(binSegmentShadowEdgeImage);

% extract edges of sample with maximum amplitude and lowest amplitude
% return is (x, ymin, ymax, 1) in homogeneous image coordinate
% this algorithm is running from 0,0 and scan start in column and row
% innner loop
sampleEdges = extractEdges(binSegmentImage);

    %Transform each half of sample edges to world unit via homography
    %transform
%     for i = 1:size(sampleEdges,1)
%         sampleHalfMinEdgesWorld(i,1:2) = homographyTransform(homographyMatrix, [sampleEdges(i,1) sampleEdges(i,2) sampleEdges(i,4)]');
%         sampleHalfMaxEdgesWorld(i,1:2) = homographyTransform(homographyMatrix, [sampleEdges(i,1) sampleEdges(i,3) sampleEdges(i,4)]');
%     end

% after getting ymax and ymin in each of x coordinate, create center plane
middlePointClound(:,1) = sampleEdges(:,1);    % x coordinate in 3D
middlePointClound(:,2) = sampleEdges(:,2);    % y min edges in 3D
middlePointClound(:,3) = sampleEdges(:,3);    % y max edge in 3D
middlePointClound(:,4) = 0;                     % z coordinate (this plane is ref plane)
% in world coordinate the image axis and world axis is converted. y=x, x=y
% of point clound
% middlePointClound(:,1) = sampleHalfMinEdgesWorld(:,2);    % x coordinate in 3D
% middlePointClound(:,2) = sampleHalfMinEdgesWorld(:,1);    % y min edges in 3D
% middlePointClound(:,3) = sampleHalfMaxEdgesWorld(:,1);    % y max edge in 3D
% middlePointClound(:,4) = 0;                     % z coordinate (this plane is ref plane)

plot3(middlePointClound(:,1), middlePointClound(:,2), middlePointClound(:,4), '.', 'LineWidth', 1);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
plot3(middlePointClound(:,1), middlePointClound(:,3), middlePointClound(:,4), '.', 'LineWidth', 1);


% find centroid point via image moment.
%which was estimated from first order image moments centroid
[centroid, headSample, tailSample] = middleLineSample(binSegmentImage);
% algorithm to find the centroid line from centroid point.
% finding reference horizontal vector in right half plane from centroid
% point.
%unitVectRefHor = calUnitVector(PD(2,:), PC(2,:));   

% finding reference horizontal vector in left half plane from centroid
% point.
n=1;     % for reference loop
    for i = int32(centroid(1)):-1:1
        centroidLine(n,1:3) = [i, centroid(2) 0];
        n = n+1;
        if(binSegmentImage(int32(centroid(2)), i) < 1)
            break;
        end
    end
    for i = int32(centroid(1)):size(binSegmentImage,2)
        centroidLine(n,1:3) = [i, centroid(2) 0];
        n = n+1;
        if(binSegmentImage(int32(centroid(2)), i) < 1)
            break;
        end
    end
% for debugging.    
%plot3([centroid(1) centroidLineR(1)], [centroid(2) centroidLineR(2)], [0 0], '*-', 'LineWidth', 1);


% for debugging.
%plot3([centroid(1) centroidLineL(1)], [centroid(2) centroidLineL(2)], [0 0], '*-', 'LineWidth', 1);
% finding centroid line by 1st order polynomial approximation
%coefficients = polyfit([double(centroidLineL(1)), double(centroidLineR(1))], [double(centroidLineL(2)), double(centroidLineR(2))], 1)
%polyMiddleLine = [sampleEdges(:,1) polyval(coefficients, sampleEdges(:,1)) zeros(size(sampleEdges(:,1),1),1)];
plot3(centroidLine(:,1), centroidLine(:,2), centroidLine(:,3), '.', 'LineWidth', 1);

% Shadow segmentation,
%---------------------------------------------------------------------
% color thresholding (HSV model)
[binSegmentShadowImage,MaskedSegmentShadowRGBImage] = shadowMask(grennScreenImage);
% filtering for only 1 largest object
binSegmentShadowImage = bwareafilt(binSegmentShadowImage, 1);
% fill hole
binSegmentShadowImage = imfill(binSegmentShadowImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentShadowImage  = medfilt2(binSegmentShadowImage ,[3 3]);
figure; imshow(binSegmentShadowImage);
