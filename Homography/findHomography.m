close all;
clear;
clc;

load('../20.09.2019/cameraParams.mat');
%imOrig = imread('../sampleImages/background/01-03-2019.JPG');
%imOrig = imread('../sampleImages/Fish Test/DSC_0332.JPG');
[filename, pathname] = uigetfile({'*.png';'*.jpg';'*.bmp'},'File Selector');
% imOrig = imread('../20.09.2019/Pillas/Picture 2019-09-20 22-07-18.PNG');
imOrig = imread(strcat(pathname, filename));

%magnification = 25;
[im, newOrigin] = undistortImage(imOrig, cameraParams, 'OutputView', 'full');
%figure; imshow(im, 'InitialMagnification', magnification);
figure; imshow(im);
title('Undistorted Image');

hold on;
colormap(gray)
axis image;
axis ij ;

dat=ginput(4);
datw = [1 1 1 1]';
Data=[dat datw]';
plot(Data, 'g+')
%w=265.6108 ; %mm 5 cm
%h=262.7235 ; %mm 5 cm
w=150.00 ; %mm 5 cm
h=150.00 ; %mm 5 cm

%paper in the real world
DataR=[ 0,0,h,h;...
       0,w,w,0;...
       1,1,1,1];
%compute the homography relating the two planes
X=diag(-DataR(1,:));
Y=diag(-DataR(2,:));
Q=Data';
A=[ Q,zeros(size(Q)),X*Q;...
        zeros(size(Q)),Q,Y*Q]
[U,S,V]=svd(A,0);
h=V(:,end);
H=reshape(h,3,3)';

% p=ginput(2)
% pw=[1 1]';
% P=[p pw]'
% 
% PR=H*P
% %Computing the distance of diameter of the coins
% delx=((PR(1,1)/PR(3,1))-(PR(1,2)/PR(3,2)));
% dely=((PR(2,1)/PR(3,1))-(PR(2,2)/PR(3,2)));
% Lenght = sqrt(delx^2+dely^2)
% plot(P,'g.-')
% text(p(1,1),p(1,2),[num2str(Lenght)]);

