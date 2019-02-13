clear;
clc;

load('cameraParams.mat');

R =[0.999754876586451,-0.002339715563935,0.022016186606441;
    0.003634721363716,0.998253387993276,-0.058965771092672;
    -0.021839769738207,0.059031339905385,0.998017196929370];

T = [-1.568660138072978e+02,-67.838402390008700,1.068218666557448e+03];

% dat = [1.592689725330621e+03,6.976576805696841e+02;
%        3.883523906408952e+03,7.096515768056961e+02;
%        3.877526958290947e+03,2.994488809766022e+03];

img = imread('undistorted.jpg');

%input from user
figure;
hold on;
colormap(gray)
imshow(img);
axis image;
axis ij ;
dat=ginput(3);
plot(dat, 'g+')


% Get the world coordinates of the corners            
worldPoints1 = pointsToWorld(cameraParams, R, T, dat);

% Compute the diameter of the coin in millimeters.
horizontalRange = worldPoints1(2, :) - worldPoints1(1, :);
verticalRange = worldPoints1(3, :) - worldPoints1(2, :);
horizontalInMillimeters = hypot(horizontalRange(1), horizontalRange(2))
verticalInMillimeters = hypot(verticalRange(1), verticalRange(2))
