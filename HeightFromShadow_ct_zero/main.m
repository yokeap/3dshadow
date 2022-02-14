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
% homographyMatrix = ...
%        [5.079053647431133e-06,-5.475219832828209e-04,0.224511299503042;
%        -5.494275672572568e-04,-1.011837140387550e-06,0.974457649498132;
%        1.174112215996749e-08,-1.921899063852394e-08,-0.005134593677710];

[filename, pathname] = uigetfile({'*.mat'},'Select session.mat','../experiment');
load(strcat(pathname, filename));



% virtual light source position
% virLightPos = [-3.313541624802504e+02,1.452013369124612e+02,4.402924880599082e+02];
% virLightPosIMG = [3.085196493165559e+03,-2.939137269817536e+03,4.361294010244374e+02];

%fish test
% virLightPos = [-3.313541624802504e+02,1.452013369124612e+02,4.402924880599082e+02];
% virLightPosIMG = [3.085196493165559e+03,-2.939137269817536e+03,4.361294010244374e+02];
 
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
% load('../calibration/cameraParams.mat');

% load background image.
%imgBg = imread('../SampleImages/background/01-03-2019.JPG');
[filename, pathname] = uigetfile({'*.png';'*.jpg';'*.bmp'},'load background image','../experiment');
imgBg = imread(strcat(pathname, filename));
% imgBg = imread('../SampleImages/Fish Test/DSC_0332.JPG');
% load sample image.
%imgSample = imread('../SampleImages/01-03-2019/DSC_0319.JPG');
[filename, pathname] = uigetfile({'*.png';'*.jpg';'*.bmp'},'load object image','../experiment');
imgBg = imread(strcat(pathname, filename));
% imgBg = imread('../SampleImages/Fish Test/DSC_0331.JPG');

%% pre-processing image procedure

% seperate object by subtracts sample image and background.
diffImage = imsubtract(rgb2gray(imgBg), rgb2gray(imgSample));
% convert to binary and filtering image.
diffImageBW = imbinarize(diffImage,0.1);
% choose largest of closed area 
diffImageBW = bwareafilt(diffImageBW, 1);
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
%[binSegmentImage, MaskedSegmentRGBImage] = createEasyMaskMango(greenScreenImage); 
% filtering for only 1 largest object
binSegmentImage = bwareafilt(binSegmentImage, 1);             
% fill hole 
binSegmentImage = imclose(binSegmentImage, strel('disk', 13));
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

% % masking sample segmented on color screen
% for i=1:size(imgSample,2)
%     for j=1:size(imgSample,1)
%         if(binSegmentImage(j,i) > 0)
%             greenScreenImage(j,i,:) = [0 0 255];
%         end      
%     end
% end

%fish testing
% masking sample segmented on color screen
for i=1:size(imgSample,2)
    for j=1:size(imgSample,1)
        if(binSegmentImage(j,i) > 0)
            diffImageBW(j,i,:) = 0;
        end      
    end
end

imshow(diffImageBW);
%% Shadow segmentation, and create middle plane of 3D

% color thresholding (HSV model)
[binSegmentShadowImage,MaskedSegmentShadowRGBImage] = createEasyMaskShadowApple2(greenScreenImage);
%[binSegmentShadowImage,MaskedSegmentShadowRGBImage] = createEasyMaskShadowMango(greenScreenImage);
% filtering for only 1 largest object
binSegmentShadowImage = bwareafilt(binSegmentShadowImage, 1);
% fill hole
binSegmentShadowImage = imclose(binSegmentShadowImage, strel('disk', 13));
binSegmentShadowImage = imfill(binSegmentShadowImage,'holes');
% apply median filtering with 10x10 kernel for smoothing edge
binSegmentShadowImage  = medfilt2(binSegmentShadowImage ,[3 3]);
figure('Name','Shadow Segmented');
imshow(binSegmentShadowImage);

% shadow edges extraction.
[sampleShadowEdgesUpper, sampleShadowEdgesLower] = extractEdges(diffImageBW);

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

%% fidning basic geometry (centroid point, central cross line)

% find centroid line by min and max of  point via image moment.
% estimating from first order image moments centroid
[centroid, headSample, tailSample] = middleLineSample(binSegmentImage);

% finding reference horizontal central line with respect with shadow ray
% direction
n=1; 
    % right half plane scan.
    for i = int32(centroid(1)):size(binSegmentImage,2)
        %centroidLineRight(n,1:3) = [i, centroid(2) 0];
        n = n+1;
        if(binSegmentImage(int32(centroid(2)), i) < 1)
            break;
        end
    end
    headSample = [i, centroid(2) 0];
    % for left half plane scan.
n=1;    
    for i = int32(centroid(1)):-1:1
        %centroidLineLeft(n,1:3) = [i, centroid(2) 0];
        n = n+1;
        if(binSegmentImage(int32(centroid(2)), i) < 1)
            break;
        end
     tailSample = [i, centroid(2) 0];
    end
    
% finding coeficient   
coefficients = polyfit([double(headSample(1)) double(tailSample(1))], [double(headSample(2)) double(tailSample(2))], 1);  
% interpolate centroid line
centroidLine = [sampleRefPointCloud(:,1) polyval(coefficients, sampleRefPointCloud(:,1)) zeros(size(sampleRefPointCloud(:,1),1),1)];
plot3(centroidLine(:,1), centroidLine(:,2), centroidLine(:,3), '.', 'LineWidth', 1);

plot3(shadowRefPointCloud(:,1), shadowRefPointCloud(:,3), shadowRefPointCloud(:,4), '.', 'LineWidth', 1);
% plot virtual light position with respect to tower c and ds
%plot3([PC(2,1) virLightPosIMG(1) PD(2,1)], [PC(2,2) virLightPosIMG(2) PD(2,2)], [PC(2,3) virLightPosIMG(3) PD(2,3)], '*-', 'LineWidth', 1, 'Color', 'r');
grid on;
hold off;

%% Transform to world unit

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

%% calculaling height from shadow

[ptHeight, sliceNumber] = calHeightFromShadow(homographyMatrix, centroidLine, binSegmentShadowImage, virLightPosIMG, virLightPos, 1000);

middleCloudLeft(:, 1) = sampleRefPointCloudWorldUpper(:, 1);
middleCloudLeft(:, 2) = sampleRefPointCloudWorldUpper(:, 2);
middleCloudLeft(:, 3) = ptHeight(:,1)/2;

middleCloudRight(:, 1) = sampleRefPointCloudWorldLower(:, 1);
middleCloudRight(:, 2) = sampleRefPointCloudWorldLower(:, 2);
middleCloudRight(:, 3) = ptHeight(:,1)/2;

% assigne height value to uppper clound
upperCloud(:, 1) = centroidLineWorld(:, 1);
upperCloud(:, 2) = centroidLineWorld(:, 2);
upperCloud(:, 3) = ptHeight(:,1) + ptHeight(:,1);

% assigne height value to lower clound
lowerCloud(:, 1) = centroidLineWorld(:, 1);
lowerCloud(:, 2) = centroidLineWorld(:, 2);
lowerCloud(:, 3) = 0;

plot3(upperCloud(:,1), upperCloud(:,2), upperCloud(:,3), '.', 'LineWidth', 1);
plot3(lowerCloud(:,1), lowerCloud(:,2), lowerCloud(:,3), '.', 'LineWidth', 1);

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
totalLength = sampleRefPointCloudWorldUpper(size(sampleRefPointCloudWorldUpper,1),2) - sampleRefPointCloudWorldUpper(1,2);
% scan along x axis 
for i=1:size(sampleRefPointCloudWorldUpper,1)
    halfSlice_L(:, 1, i) = [upperCloud(i, 1), middleCloudLeft(i, 1), lowerCloud(i, 1)];
    halfSlice_L(:, 2, i) = [upperCloud(i, 3), middleCloudLeft(i, 3), lowerCloud(i, 3)];
    
    halfSlice_R(:, 1, i) = [lowerCloud(i, 1), middleCloudRight(i, 1), upperCloud(i, 1)];
    halfSlice_R(:, 2, i) = [lowerCloud(i, 3), middleCloudRight(i, 3), upperCloud(i, 3)];
    
    
     %spline interpolation
     sliceModel_L(:,2:3,i) = splineEstimate(halfSlice_L(:,:,i));
     sliceModel_R(:,2:3,i) = splineEstimate(halfSlice_R(:,:,i));
     
     sliceModel(:,2:3,i) = [sliceModel_L(:,2:3,i); sliceModel_R(:,2:3,i)];
    
     
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
     sliceModel(:,1,i)= sampleRefPointCloudWorldUpper(i,2);
     sliceModelArea(i,1) = sampleRefPointCloudWorldUpper(i,2);
     
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

figure('Name','3D Point Cloud');
%plot3(leftSliceX,leftSliceY,leftSliceZ, 'r.');
plot3(sliceX,sliceY,sliceZ, 'r.');
hold on
%plot3(rightSliceX,rightSliceY,rightSliceZ, 'b.');
xlabel('X (mm)','FontSize',20,'FontWeight','bold');
ylabel('Y (mm)','FontSize',20,'FontWeight','bold');
zlabel('Height (mm)','FontSize',20,'FontWeight','bold');
title(['Volume = ' num2str(totalVolume/1000) ' cm^3'],'FontSize',20,'FontWeight','bold','Color','b');
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

% figure('Name', '2D Cross section test');
% plot(sliceModelRaw(:,2,sliceNumber), sliceModelRaw(:,3,sliceNumber), 'r.');
