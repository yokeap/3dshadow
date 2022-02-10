close all;
clear;
clc;
%%  - load calibrated variable
%   - load raw image and background image. 

% homographyMatrix = ...
%     [-1.923406815982628e-06,-5.164347311016744e-04,0.192943934722844;
%     -5.617555953485216e-04,-1.940433835436591e-06,0.981195951656823;
%     -7.321030736267486e-09,-1.093315495754124e-08,-0.005153661250271];

%fish test
homographyMatrix = ...
       [5.079053647431133e-06,-5.475219832828209e-04,0.224511299503042;
       -5.494275672572568e-04,-1.011837140387550e-06,0.974457649498132;
       1.174112215996749e-08,-1.921899063852394e-08,-0.005134593677710];

% virtual light source position
% virLightPos = [-3.313541624802504e+02,1.452013369124612e+02,4.402924880599082e+02];
% virLightPosIMG = [3.085196493165559e+03,-2.939137269817536e+03,4.361294010244374e+02];

%fish test
virLightPos = [-3.382978646119368e+02,1.228094975199053e+02,2.164641969419283e+02];
virLightPosIMG = [2.907071993662571e+03,-2.682554285912778e+03,2.164641969419283e+02];
 
% azimuth triangle (XY) or top view in radian
% azimuth_a = 1.3456;
% azimuth_b = 1.2907;
% azimuth_c = 0.5053;     % this angle is used for ray casting in top view
 
% elevation triangle (XZ) or side view in radian
% elevation_a = 0.3231;
% elevation_b = 2.8158;
% elevation_c = 0.0027;

% elevation_a = 0.8267;
% elevation_b = 1.6186;
% elevation_c = 0.6963;

 % load intrinsic camera calibrated parameters
%load('../calibration/cameraParams.mat');

% load background image.
%imgBg = imread('../SampleImages/background/01-03-2019.JPG');
imgBg = imread('../SampleImages/Fish Test/DSC_0332.JPG');
% load sample image.
%imgSample = imread('../SampleImages/01-03-2019/DSC_0319.JPG');
imgSample = imread('../SampleImages/Fish Test/DSC_0331.JPG');

%% pre-processing image procedure

% seperate object by subtracts sample image and background.
diffImage = imsubtract(rgb2gray(imgBg), rgb2gray(imgSample));
% convert to binary and filtering image.
diffImageBW = imbinarize(diffImage,0.1);
% choose largest of calosed area 
diffImageBW = bwareafilt(diffImageBW, 1);
figure('Name','Original');
imshow(diffImageBW);
% fill hole of area
diffImageBW = imfill(diffImageBW,'holes');
figure('Name','Differential Image');
imshow(diffImageBW);

% create blank blackgroumd image (background color is used to seperated the object).
% ie. apple used green screen, mango used red screen.
greenScreenImage=zeros(size(imgSample,1), size(imgSample,2), 3,'uint8');
greenScreenImage(:,:,3) = 255;         % for green (channel 2)
%greenScreenImage(:,:,2) = 255;         % for green (channel 2)
%greenScreenImage(:,:,1) = 255;          % for red (channel 1)

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

%% Sample segmentation, and create middle plane of 3D

% color thresholding (HSV model)
[binSegmentImage, MaskedSegmentRGBImage] = createEasyMaskFish(greenScreenImage);     
figure('Name','MaskedSegmentRGBImage');
imshow(MaskedSegmentRGBImage);
%[binSegmentImage, MaskedSegmentRGBImage] = createEasyMaskMango(greenScreenImage); 
% filtering for only 1 largest object
binSegmentImage = bwareafilt(binSegmentImage, 1);             
% fill hole 
%binSegmentImage = imclose(binSegmentImage, strel('disk', 13));
binSegmentImage = imclose(binSegmentImage, strel('disk', 300));
binSegmentImage = imfill(binSegmentImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentImage  = medfilt2(binSegmentImage ,[15 15]);
figure('Name','Sample Segmented');
imshow(binSegmentImage);

%% Shadow segmentation, and create middle plane of 3D

%fish testing
% masking sample segmented on color screen
for i=1:size(imgSample,2)
    for j=1:size(imgSample,1)
        if(binSegmentImage(j,i) > 0)
            diffImageBW(j,i,:) = 0;
        end      
    end
end
% filtering for only 1 largest object
binSegmentShadowImage = bwareafilt(diffImageBW, 2);
% fill hole
binSegmentShadowImage = imclose(binSegmentShadowImage, strel('disk', 13));
binSegmentShadowImage = imfill(binSegmentShadowImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentShadowImage  = medfilt2(binSegmentShadowImage ,[3 3]);
figure('Name','Shadow Segmented');
imshow(binSegmentShadowImage);

%% calculaling height from shadow

[middleCloudLeft, middleCloudRight, shadowHeight] = reconstruct(homographyMatrix, binSegmentImage, MaskedSegmentRGBImage, binSegmentShadowImage, virLightPosIMG, virLightPos, 2000);

maxVal = max(shadowHeight(:, 3));
% middleCloudLeft(:, 1) = sampleRefPointCloudWorldUpper(:, 1);
% middleCloudLeft(:, 2) = sampleRefPointCloudWorldUpper(:, 2);
middleCloudLeft(:, 3) = 0;

% middleCloudRight(:, 1) = sampleRefPointCloudWorldLower(:, 1);
% middleCloudRight(:, 2) = sampleRefPointCloudWorldLower(:, 2);
middleCloudRight(:, 3) = 0;

% assigne height value to uppper clound
upperCloud(:, 1) = shadowHeight(:, 1);
upperCloud(:, 2) = shadowHeight(:, 2);
upperCloud(:, 3) = shadowHeight(:, 3) / 2;

% assigne height value to lower clound
lowerCloud(:, 1) = shadowHeight(:, 1);
lowerCloud(:, 2) = shadowHeight(:, 2);
lowerCloud(:, 3) = upperCloud(:, 3) * -1;
%lowerCloud(:, 3) = maxVal - skeleton(:, 3);



figure('Name','3D Points Cloud (World)');
hold on
plot3(middleCloudLeft(:,1), middleCloudLeft(:,2), middleCloudLeft(:,3), '.', 'LineWidth', 1);
plot3(middleCloudRight(:,1), middleCloudRight(:,2), middleCloudRight(:,3), '.', 'LineWidth', 1);
plot3(upperCloud(:,1), upperCloud(:,2), upperCloud(:,3), '.', 'LineWidth', 1);
plot3(lowerCloud(:,1), lowerCloud(:,2), lowerCloud(:,3), '.', 'LineWidth', 1);
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on;


%% volume estimating section

% allocate the variable
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
totalLength = middleCloudLeft(size(middleCloudLeft,1),2) - middleCloudLeft(1,2);
% scan along x axis 
for i=1:size(middleCloudLeft,1)
    
%     halfSlice_R(:, 1, i) = [lowerCloud(i, 1), middleCloudRight(i, 1), upperCloud(i, 1)];
%     halfSlice_R(:, 2, i) = [lowerCloud(i, 3), middleCloudRight(i, 3), upperCloud(i, 3)];
    
    if abs(middleCloudRight(i, 1) - upperCloud(i, 1)) < abs(middleCloudLeft(i, 1) - upperCloud(i, 1))
%         mirrorUpperXPosFrom_R = middleCloudLeft(i, 1) + abs(middleCloudRight(i, 1) - upperCloud(i, 1));
        halfSliceUpper(:, 1, i) = [middleCloudRight(i, 1), upperCloud(i, 1), mirrorUpperXPosFrom_R, middleCloudLeft(i, 1)];
        halfSliceUpper(:, 2, i) = [middleCloudRight(i, 3), upperCloud(i, 3), upperCloud(i, 3), middleCloudLeft(i, 3)];
    else
        mirrorUpperXPosFrom_L = middleCloudRight(i, 1) - abs(middleCloudLeft(i, 1) - upperCloud(i, 1));
        halfSliceUpper(:, 1, i) = [middleCloudRight(i, 1), mirrorUpperXPosFrom_L, upperCloud(i, 1), middleCloudLeft(i, 1)];
        halfSliceUpper(:, 2, i) = [middleCloudRight(i, 3), upperCloud(i, 3), upperCloud(i, 3), middleCloudLeft(i, 3)];
    end
    
    
    if abs(middleCloudRight(i, 1) - lowerCloud(i, 1)) < abs(middleCloudLeft(i, 1) - lowerCloud(i, 1))
        mirrorLowerXPosFrom_R = middleCloudLeft(i, 1) + abs(middleCloudRight(i, 1) - lowerCloud(i, 1));
        halfSliceLower(:, 1, i) = [middleCloudLeft(i, 1), mirrorLowerXPosFrom_R, lowerCloud(i, 1), middleCloudRight(i, 1)];
        halfSliceLower(:, 2, i) = [middleCloudLeft(i, 3), lowerCloud(i, 3), lowerCloud(i, 3), middleCloudRight(i, 3)];
    else
        mirrorLowerXPosFrom_L = middleCloudRight(i, 1) - abs(middleCloudLeft(i, 1) - lowerCloud(i, 1));
        halfSliceLower(:, 1, i) = [middleCloudLeft(i, 1), lowerCloud(i, 1), mirrorLowerXPosFrom_L, middleCloudRight(i, 1)];
        halfSliceLower(:, 2, i) = [middleCloudLeft(i, 3), lowerCloud(i, 3), lowerCloud(i, 3), middleCloudRight(i, 3)];
    end
    
%     halfSliceUpper(:, 1, i) = [middleCloudRight(i, 1), upperCloud(i, 1), mirrorUpperXPosFrom_R, middleCloudLeft(i, 1)];
%     halfSliceUpper(:, 2, i) = [middleCloudRight(i, 3), upperCloud(i, 3), upperCloud(i, 3), middleCloudLeft(i, 3)];
%     
%     halfSliceLower(:, 1, i) = [middleCloudLeft(i, 1), mirrorLowerXPosFrom_R, lowerCloud(i, 1), middleCloudRight(i, 1)];
%     halfSliceLower(:, 2, i) = [middleCloudLeft(i, 3), lowerCloud(i, 3), lowerCloud(i, 3), middleCloudRight(i, 3)];
    
    
    
    
%     mirrorUpperCloudFrom_R = [mirrorXPosFrom_R, upperCloud(i, 3)]; 
%     mirrorLowerCloudFrom_R = [mirrorXPosFrom_R, lowerCloud(i, 3)]; 
    
%     halfSlice_L(:, 1, i) = [mirrorUpperXPosFrom_R, middleCloudLeft(i, 1), mirrorLowerXPosFrom_R];
%     halfSlice_L(:, 2, i) = [upperCloud(i, 3), middleCloudLeft(i, 3), lowerCloud(i, 3)];
    
     %spline interpolation
     sliceModel_Upper(:,2:3,i) = splineEstimate(halfSliceUpper(:,:,i));
     sliceModel_Lower(:,2:3,i) = splineEstimate(halfSliceLower(:,:,i));
     
     sliceModel(:,2:3,i) = [sliceModel_Upper(:,2:3,i); sliceModel_Lower(:,2:3,i)];
    
     
%      sliceModelRaw(:,2:3,i) = [upperCloud(i, 1), upperCloud(i, 3);...
%                           sampleRefPointCloudWorldUpper(i, 1), sampleRefPointCloudWorldUpper(i, 3);...
%                           lowerCloud(i, 1), lowerCloud(i, 3);...
%                            %lowerCloud(i, 1), lowerCloud(i, 3);...
%                           sampleRefPointCloudWorldLower(i, 1), sampleRefPointCloudWorldLower(i, 3);...
%                           upperCloud(i, 1), upperCloud(i, 3)];
                      
%      sliceModelRaw(:,2:3,i) = [sampleRefPointCloudWorldUpper(i, 1), ptHeight(i);...
%                                sampleRefPointCloudWorldLower(i, 1), ptHeight(i);...
%                                sampleRefPointCloudWorldLower(i, 1), ptHeight(i)*-1;...
%                                sampleRefPointCloudWorldUpper(i, 1), ptHeight(i)*-1;...
%                                sampleRefPointCloudWorldUpper(i, 1), ptHeight(i)]; %closed point
     

                           
      %sliceModel(:,2:3,i) = splineEstimate(sliceModelRaw(:,2:3,i));                
     sliceModelArea(i,2) = polyarea(sliceModel(:,2,i),sliceModel(:,3,i));
        
     % append y (Depth) coordinate to column 1
     %sliceModel_L(:,1,i) = sampleRefPointCloudWorldUpper(i,2);
     %sliceModel_R(:,1,i) = sampleRefPointCloudWorldUpper(i,2);
     %sliceModelRaw(:,1,i)= sampleRefPointCloudWorldUpper(i,2);
     sliceModel(:,1,i)= middleCloudLeft(i,2);
     sliceModelArea(i,1) = middleCloudLeft(i,2);
     
     if i>1
%          sumArea = sliceModelArea(i,2) + sliceModelArea(i-1,2);
           sumVolume = sliceModelArea(i,2) * (sliceModelArea(i,1) - sliceModelArea(i-1,1));
%          sumVolume = sliceModelArea(i,2) * 0.1165;
           totalVolume = totalVolume + sumVolume;

%           sA = ((sliceModelArea(i,2) + sliceModelArea(i-1,2))/2) * (sliceModelArea(i,1) - sliceModelArea(i-1,1)) + ...
%                  (sliceModelArea(1,2)/2)*totalLength + (sliceModelArea(size(sliceModelArea,1),2)/2)*totalLength;
%            totalVolume = totalVolume + sA;
     end

    sliceX = [sliceX; sliceModel(:,1,i)];
    sliceY = [sliceY; sliceModel(:,2,i)]; 
    sliceZ = [sliceZ; sliceModel(:,3,i)];
end

figure('Name', 'Upper Slice');
hold on;
plot(halfSliceUpper(:, 1, 50), halfSliceUpper(:, 2, 50), 'r.');
grid on;

% figure('Name', 'Right Slice');
% hold on;
% plot(halfSliceUpper(:, 1, 50), halfSliceUpper(:, 2, 50), 'r.');
% grid on;
% 
% figure('Name', 'Left Slice');
% hold on;
% plot(halfSliceUpper(:, 1, 50), halfSliceLower(:, 2, 50), 'r.');
% grid on;

figure('Name','3D Point Cloud');
%plot3(leftSliceX,leftSliceY,leftSliceZ, 'r.');
plot3(sliceX,sliceY,sliceZ);
hold on;

%plot3(rightSliceX,rightSliceY,rightSliceZ, 'b.');
xlabel('X (mm)','FontSize',20,'FontWeight','bold');
ylabel('Y (mm)','FontSize',20,'FontWeight','bold');
zlabel('Height (mm)','FontSize',20,'FontWeight','bold');
title(['Volume = ' num2str(totalVolume/1000) ' cm^3'],'FontSize',20,'FontWeight','bold','Color','b');
grid on;

c = 1:numel(sliceZ);      %# colors
h = surface([sliceX(:), sliceX(:)], [sliceY(:), sliceY(:)], [sliceZ(:), sliceZ(:)], ...
    [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
colormap( jet(numel(sliceZ)) )

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
grid on
% figure('Name', '2D Cross section test');
% plot(sliceModelRaw(:,2,sliceNumber), sliceModelRaw(:,3,sliceNumber), 'r.');
