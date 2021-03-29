close all;
clear;
clc;

homographyMatrix = ...
[ -0.000001171027520  -0.000569119373972   0.402795558322037
  -0.000571860777008  -0.000007045512653   0.915276425974051
  -0.000000000072476  -0.000000017684098  -0.004914395139275];

% virtual light source position
virLightPos = [-2.298030767150990e+02,1.469862104150302e+02,4.943304818724096e+02];
virLightPosIMG = [2.868181414399634e+03,-1.279352211841624e+03,4.871044729212234e+02];

% Reference tower position in image coordinate.
PA = [1.591619136460555e+03,1.426180970149254e+03;7.014912046908303e+02,9.517694562899778e+02;1,1]';
  
PB = [3.873817430703625e+03,4.005319562899787e+03;7.014912046908303e+02,9.517694562899778e+02;1,1]';

PC = [3.873817430703625e+03,3.992593550106611e+03;3.174579690831556e+03,3.700588219616204e+03;1,1]';

PD = [1.570409115138593e+03,1.413454957356077e+03;3.174579690831556e+03,3.713314232409381e+03;1,1]';

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
% elevation_a = 0.3231;
% elevation_b = 2.8158;
% elevation_c = 0.0027;

elevation_a = 0.8267;
elevation_b = 1.6186;
elevation_c = 0.6963;

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
figure('Name','Differential Image');
imshow(diffImageBW);

% create blank blackgroumd image (background color is according to color of object).
greenScreenImage=zeros(size(imgSample,1), size(imgSample,2), 3,'uint8');
greenScreenImage(:,:,2) = 255;

% masking the input image from masked image that come from subtracted image.
% by checking boolean of image.
% *note: matlab direct accessing pixel is (y,x), ie img(y,x).
for i=1:size(imgSample,2)
    for j=1:size(imgSample,1)
        if(diffImageBW(j,i) > 0)
            greenScreenImage(j,i,:) = imgSample(j,i,:);
        end      
    end
end

figure('Name','Green Screen Image');
imshow(greenScreenImage);

% Sample segmentation,
%---------------------------------------------------------------------
% color thresholding (HSV model)
[binSegmentImage, MaskedSegmentRGBImage] = appleMask(greenScreenImage);        
% filtering for only 1 largest object
binSegmentImage = bwareafilt(binSegmentImage, 1);             
% fill hole
binSegmentImage = imfill(binSegmentImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentImage  = medfilt2(binSegmentImage ,[3 3]);
figure('Name','Sample Segmented');
imshow(binSegmentImage);

% extract uppper and lower edges of sample and return
% in homogeneous image coordinate
% this algorithm is running from 0,0 and scan start in column and row
[sampleEdgesUpper, sampleEdgesLower] = extractEdges(binSegmentImage);

% create center plane point Cloud.
sampleRefPointCloud(:,1) = sampleEdgesUpper(:,1);    % x coordinate in 3D
sampleRefPointCloud(:,2) = sampleEdgesUpper(:,2);    % y min edges in 3D
sampleRefPointCloud(:,3) = sampleEdgesLower(:,2);    % y max edge in 3D
sampleRefPointCloud(:,4) = 0;                     % z coordinate (this plane is ref plane)

% Shadow segmentation,
%---------------------------------------------------------------------
% color thresholding (HSV model)
[binSegmentShadowImage,MaskedSegmentShadowRGBImage] = shadowMask(greenScreenImage);
% filtering for only 1 largest object
binSegmentShadowImage = bwareafilt(binSegmentShadowImage, 1);
% fill hole
binSegmentShadowImage = imfill(binSegmentShadowImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentShadowImage  = medfilt2(binSegmentShadowImage ,[3 3]);
figure('Name','Shadow Segmented');
imshow(binSegmentShadowImage);

% shadow edges extraction.
[sampleShadowEdgesUpper, sampleShadowEdgesLower] = extractEdges(binSegmentShadowImage);

shadowRefPointCloud(:,1) = sampleShadowEdgesUpper(:,1);    % x coordinate in 3D
shadowRefPointCloud(:,2) = sampleShadowEdgesUpper(:,2);    % y min edges in 3D
shadowRefPointCloud(:,3) = sampleShadowEdgesLower(:,2);    % y max edge in 3D
shadowRefPointCloud(:,4) = 0;                         % z coordinate (this plane is ref plane)

figure('Name','3D Points Cloud');
plot3(sampleRefPointCloud(:,1), sampleRefPointCloud(:,2), sampleRefPointCloud(:,4), '.', 'LineWidth', 1);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
plot3(sampleRefPointCloud(:,1), sampleRefPointCloud(:,3), sampleRefPointCloud(:,4), '.', 'LineWidth', 1);

% find centroid line by min and max of  point via image moment.
% estimating from first order image moments centroid
[centroid, headSample, tailSample] = middleLineSample(binSegmentImage);

% finding reference horizontal line centroid point.
n=1;     % for reference loop
    %right half plane.
    for i = int32(centroid(1)):size(binSegmentImage,2)
        %centroidLineRight(n,1:3) = [i, centroid(2) 0];
        n = n+1;
        if(binSegmentImage(int32(centroid(2)), i) < 1)
            break;
        end
    end
    headSample = [i, centroid(2) 0];
    % for left half plane.
n=1;    
    for i = int32(centroid(1)):-1:1
        %centroidLineLeft(n,1:3) = [i, centroid(2) 0];
        n = n+1;
        if(binSegmentImage(int32(centroid(2)), i) < 1)
            break;
        end
     tailSample = [i, centroid(2) 0];
    end
    
coefficients = polyfit([double(headSample(1)) double(tailSample(1))], [double(headSample(2)) double(tailSample(2))], 1);    
centroidLine = [sampleRefPointCloud(:,1) polyval(coefficients, sampleRefPointCloud(:,1)) zeros(size(sampleRefPointCloud(:,1),1),1)];
plot3(centroidLine(:,1), centroidLine(:,2), centroidLine(:,3), '.', 'LineWidth', 1);

plot3(shadowRefPointCloud(:,1), shadowRefPointCloud(:,3), shadowRefPointCloud(:,4), '.', 'LineWidth', 1);
% plot virtual light position with respect to tower c and ds
%plot3([PC(2,1) virLightPosIMG(1) PD(2,1)], [PC(2,2) virLightPosIMG(2) PD(2,2)], [PC(2,3) virLightPosIMG(3) PD(2,3)], '*-', 'LineWidth', 1, 'Color', 'r');
grid on;
hold off;

% transform image coordinate to world unit
for i = 1:size(sampleRefPointCloud, 1)
    sampleRefPointCloudWorldUpper(i,:) = [homographyTransform(homographyMatrix, [sampleEdgesUpper(i, 1:2) 1]'); 0];
    sampleRefPointCloudWorldLower(i,:) = [homographyTransform(homographyMatrix, [sampleEdgesLower(i, 1:2) 1]'); 0];
    centroidLineWorld(i,:) = [homographyTransform(homographyMatrix, [centroidLine(i, 1:2) 1]'); 0];
end

for i = 1:size(shadowRefPointCloud, 1)
    shadowRefPointCloudWorldLower(i,:) = [homographyTransform(homographyMatrix, [sampleShadowEdgesLower(i, 1:2) 1]'); 0];
end

figure('Name','3D Points Cloud (World Unit)');
plot3(sampleRefPointCloudWorldUpper(:,1), sampleRefPointCloudWorldUpper(:,2), sampleRefPointCloudWorldUpper(:,3), '.', 'LineWidth', 1);
hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
plot3(sampleRefPointCloudWorldLower(:,1), sampleRefPointCloudWorldLower(:,2), sampleRefPointCloudWorldLower(:,3), '.', 'LineWidth', 1);
%plot3(shadowRefPointCloudWorldLower(:,1), shadowRefPointCloudWorldLower(:,2), shadowRefPointCloudWorldLower(:,3), '.', 'LineWidth', 1);

[ptHeight, sliceNumber] = calHeightFromShadow(homographyMatrix, centroidLine, binSegmentShadowImage, virLightPosIMG, virLightPos, 1000);

% middleCloudLower(:, 1) = sampleRefPointCloudWorldLower(:, 1);
% middleCloudLower(:, 2) = sampleRefPointCloudWorldLower(:, 2);
% middleCloudLower(:, 3) = sampleRefPointCloudWorldLower(:, 3);

upperCloud(:, 1) = centroidLineWorld(:, 1);
upperCloud(:, 2) = centroidLineWorld(:, 2);
upperCloud(:, 3) = ptHeight(:,1);

lowerCloud(:, 1) = centroidLineWorld(:, 1);
lowerCloud(:, 2) = centroidLineWorld(:, 2);
lowerCloud(:, 3) = ptHeight(:,1) * -1;

plot3(upperCloud(:,1), upperCloud(:,2), upperCloud(:,3), '.', 'LineWidth', 1);
plot3(lowerCloud(:,1), lowerCloud(:,2), lowerCloud(:,3), '.', 'LineWidth', 1);

sliceX = [];
sliceY = [];
sliceZ = [];
sliceModel_L = [];
sliceModel_R = [];
leftSliceX = [];
leftSliceY = [];
leftSliceZ = [];
rightSliceX = [];
rightSliceY = [];
rightSliceZ = [];
sumArea = 0;
sumVolume = 0;
totalVolume = 0;
totalLength = sampleRefPointCloudWorldUpper(size(sampleRefPointCloudWorldUpper,1),2) - sampleRefPointCloudWorldUpper(1,2);
for i=1:size(sampleRefPointCloudWorldUpper,1)
    halfSlice_L(:, 1, i) = [upperCloud(i, 1), sampleRefPointCloudWorldUpper(i, 1), lowerCloud(i, 1)];
    halfSlice_L(:, 2, i) = [upperCloud(i, 3), sampleRefPointCloudWorldUpper(i, 3), lowerCloud(i, 3)];
    
    halfSlice_R(:, 1, i) = [lowerCloud(i, 1), sampleRefPointCloudWorldLower(i, 1), upperCloud(i, 1)];
    halfSlice_R(:, 2, i) = [lowerCloud(i, 3), sampleRefPointCloudWorldLower(i, 3), upperCloud(i, 3)];
    
%     halfSlice_L(:, 1, i) = [sampleRefPointCloudWorldUpper(i, 1), upperCloud(i, 1), sampleRefPointCloudWorldLower(i, 1)];
%     halfSlice_L(:, 2, i) = [sampleRefPointCloudWorldUpper(i, 3), upperCloud(i, 3), sampleRefPointCloudWorldLower(i, 3)];
%     
%     halfSlice_R(:, 1, i) = [sampleRefPointCloudWorldUpper(i, 1), lowerCloud(i, 1), sampleRefPointCloudWorldLower(i, 1)];
%     halfSlice_R(:, 2, i) = [sampleRefPointCloudWorldUpper(i, 3), lowerCloud(i, 3), sampleRefPointCloudWorldLower(i, 3)];
    

     % spline interpolation
%      sliceModel_L(:,2:3,i) = splineEstimate(halfSlice_L(:,:,i));
%      sliceModel_R(:,2:3,i) = splineEstimate(halfSlice_R(:,:,i));
     
     %sliceModel(:,:,i) = [sliceModel_L(:,2:3,i); sliceModel_R(:,2:3,i)];
     sliceModelRaw(:,2:3,i) = [upperCloud(i, 1), upperCloud(i, 3);...
                          sampleRefPointCloudWorldUpper(i, 1), sampleRefPointCloudWorldUpper(i, 3);...
                          lowerCloud(i, 1), lowerCloud(i, 3);...
                           %lowerCloud(i, 1), lowerCloud(i, 3);...
                          sampleRefPointCloudWorldLower(i, 1), sampleRefPointCloudWorldLower(i, 3);...
                          upperCloud(i, 1), upperCloud(i, 3)];
                      
     sliceModel(:,2:3,i) = splineEstimate(sliceModelRaw(:,2:3,i));                
     sliceModelArea(i,2) = polyarea(sliceModel(:,2,i),sliceModel(:,3,i));
        
     % append y (Depth) coordinate to column 1
%      sliceModel_L(:,1,i) = sampleRefPointCloudWorldUpper(i,2);
%      sliceModel_R(:,1,i) = sampleRefPointCloudWorldUpper(i,2);
     sliceModelRaw(:,1,i)= sampleRefPointCloudWorldUpper(i,2);
     sliceModel(:,1,i)= sampleRefPointCloudWorldUpper(i,2);
     sliceModelArea(i,1) = sampleRefPointCloudWorldUpper(i,2);
     
      %sumArea = sumArea + sliceModelArea(i,2);
      %sumH = sumH + sliceModelArea(i,1);
%      sumVolume = sliceModelArea(i,2) * 0.1165;
%       sumVolume = sliceModelArea(i,2) * 1.165;
%       totalVolume = totalVolume + sumVolume;
     if i>1
%          sumArea = sliceModelArea(i,2) + sliceModelArea(i-1,2);
           sumVolume = sliceModelArea(i,2) * (sliceModelArea(i,1) - sliceModelArea(i-1,1));
%          sumVolume = sliceModelArea(i,2) * 0.1165;
           totalVolume = totalVolume + sumVolume;

%           sA = ((sliceModelArea(i,2) + sliceModelArea(i-1,2))/2) * (sliceModelArea(i,1) - sliceModelArea(i-1,1)) + ...
%                  (sliceModelArea(1,2)/2)*totalLength + (sliceModelArea(size(sliceModelArea,1),2)/2)*totalLength;
%            totalVolume = totalVolume + sA;
     end
%      
    %scatter3(sliceModel_R(:,1,i), sliceModel_R(:,2,i), sliceModel_R(:,3,i), 2)
%     leftSliceX = [leftSliceX; sliceModel_L(:,1,i)];
%     leftSliceY = [leftSliceY; sliceModel_L(:,2,i)];
%     leftSliceZ = [leftSliceZ; sliceModel_L(:,3,i)];
%     
%     %scatter3(sliceModel_R(:,1,i), sliceModel_R(:,2,i), sliceModel_R(:,3,i), 2)
%     rightSliceX = [rightSliceX; sliceModel_R(:,1,i)];
%     rightSliceY = [rightSliceY; sliceModel_R(:,2,i)];
%     rightSliceZ = [rightSliceZ; sliceModel_R(:,3,i)];

    sliceX = [sliceX; sliceModel(:,1,i)];
    sliceY = [sliceY; sliceModel(:,2,i)]; 
    sliceZ = [sliceZ; sliceModel(:,3,i)];
end
% sumVolume = sumArea * sumH;


figure('Name','3D Point Cloud');
%plot3(leftSliceX,leftSliceY,leftSliceZ, 'r.');
plot3(sliceX,sliceY,sliceZ, 'r.');
hold on
%plot3(rightSliceX,rightSliceY,rightSliceZ, 'b.');
xlabel('X (mm)','FontSize',20,'FontWeight','bold');
ylabel('Y (mm)','FontSize',20,'FontWeight','bold');
zlabel('Height (mm)','FontSize',20,'FontWeight','bold');
title(['Volume = ' num2str(totalVolume) ' mm^3'],'FontSize',20,'FontWeight','bold','Color','b');
grid on;

sliceNumber = 50;

% figure('Name','2D Half Left Slicing');
% plot([sliceModel_L(1,1,sliceNumber) sliceModel_L(2,1,sliceNumber) sliceModel_L(3,1,sliceNumber)], [sliceModel_L(1,2,sliceNumber) sliceModel_L(2,2,sliceNumber) sliceModel_L(3,2,sliceNumber)], 'o-', 'LineWidth', 1);
% hold on;
% plot([sliceModel_R(1,1,sliceNumber) sliceModel_R(2,1,sliceNumber) sliceModel_R(3,1,sliceNumber)], [sliceModel_R(1,2,sliceNumber) sliceModel_R(2,2,sliceNumber) sliceModel_R(3,2,sliceNumber)], '*-', 'LineWidth', 1);


% figure('Name','2D Half Left Slicing');
% plot(sliceModel_L(:,2,sliceNumber), sliceModel_L(:,3,sliceNumber), 'r.');
% hold on;
% plot(sliceModel_R(:,2,sliceNumber), sliceModel_R(:,3,sliceNumber), 'b.');

figure('Name', '2D Cross section');
title(['Area = ' num2str(sliceModelArea(sliceNumber,2))]);
hold on;
plot(sliceModel(:,2,sliceNumber), sliceModel(:,3,sliceNumber), 'r.');
