
% 3D images from shadow (shape from shadow).
% Beta version 1.
% Author, Siwakorn Sukprasertchai
% *********************************************
% short note:
% * using calibrated intrinsic parameters to get undistorted the image via
%   radial and tangential model. Using computed homography matrix to
%   estimates interesting point and transforming to the world coordinate.
% -----------------------------------------------------------------------

clear;
close all;
clc;

% homographyMatrix = ...
% [ -0.000001171027520  -0.000569119373972   0.402795558322037
%   -0.000571860777008  -0.000007045512653   0.915276425974051
%   -0.000000000072476  -0.000000017684098  -0.004914395139275];

% homographyMatrix = ...
%     [-1.923406815982628e-06,-5.164347311016744e-04,0.192943934722844;
%     -5.617555953485216e-04,-1.940433835436591e-06,0.981195951656823;
%     -7.321030736267486e-09,-1.093315495754124e-08,-0.005153661250271];

% % load instrincsic camera calibrated parameters
% load('../calibration/cameraParams.mat');
% % load computed homography matrix

% load intrinsic, extrinsic and homography matrix
load('../20.09.2019/homographyParams.mat');

homographyMatrix = H;

% load distorted image
%imOrig = imread('../SampleImages/background/01-03-2019.JPG');
imOrig = imread('../20.09.2019/Pillas/Picture 2019-09-20 22-07-18.PNG');
% get undistorting image and offset origin
[im, newOrigin] = undistortImage(imOrig, cameraParams, 'OutputView', 'full');
figure; imshow(im); title('Undistorted Image');
axis image;
axis ij ;

% get the coordinate of shadow ray
hold on;
input_A = [ginput(2)'; [1 1]];
input_B = [ginput(2)'; [1 1]]; 
input_C = [ginput(2)'; [1 1]];
input_D = [ginput(2)'; [1 1]];

% input_A = [1.742544871794872e+03,1.560673076923077e+03;3.653525641025626e+02,8.222499999999986e+02];
% input_B = [4.200032051282052e+03,4.350852564102564e+03;3.653525641025626e+02,8.222499999999986e+02];
% input_C = [4.200032051282052e+03,4.350852564102564e+03;3.004711538461538e+03,3.798737179487179e+03];
% input_D = [1.742544871794872e+03,1.560673076923077e+03;3.022455128205128e+03,3.820916666666666e+03];

% input_B = nearestPointY(input_A,input_B);
% input_C = nearestPointX(input_B,input_C);
% input_D = nearestPointX(input_A,input_D);
% input_D = nearestPointY(input_C,input_D);

% input_A = [input_A; [1 1]];
% input_B = [input_B; [1 1]];
% input_C = [input_C; [1 1]];
% input_D = [input_D; [1 1]];

plot(input_A(1,1:2), input_A(2, 1:2), '*-');
hold on;
plot(input_B(1,1:2), input_B(2, 1:2), '*-');
plot(input_C(1,1:2), input_C(2, 1:2), '*-');
plot(input_D(1,1:2), input_D(2, 1:2), '*-');
hold off;

P1=homographyTransform(homographyMatrix,input_A(:,1));
S1=homographyTransform(homographyMatrix,input_A(:,2));

P2=homographyTransform(homographyMatrix,input_B(:,1));
S2=homographyTransform(homographyMatrix,input_B(:,2));

P3=homographyTransform(homographyMatrix,input_C(:,1));
S3=homographyTransform(homographyMatrix,input_C(:,2));

P4=homographyTransform(homographyMatrix,input_D(:,1));
S4=homographyTransform(homographyMatrix,input_D(:,2));

figure;
plot([P1(1,1) S1(1,1)],[P1(2,1) S1(2,1)],'*-', 'Color', 'r');
hold on;
plot([P2(1,1) S2(1,1)],[P2(2,1) S2(2,1)],'*-', 'Color', 'g');
plot([P3(1,1) S3(1,1)],[P3(2,1) S3(2,1)],'*-', 'Color', 'b');
plot([P4(1,1) S4(1,1)],[P4(2,1) S4(2,1)],'*-', 'Color', 'y');
hold off;

% rearrange the matrix
A = [P1 S1 P1; [0 0 52]]';
B = [P2 S2 P2; [0 0 52]]';
C = [P3 S3 P3; [0 0 52]]';
D = [P4 S4 P4; [0 0 52]]';


% plot vector from shadow edges to top of object
figure('Name','Inverese Ray Tracing');
hold on;
plot3(A(2:3,1), A(2:3,2), A(2:3,3), 'o-','LineWidth',5, 'Color', 'r');
plot3(B(2:3,1), B(2:3,2), B(2:3,3),'o-','LineWidth',5, 'Color', 'g');
plot3(C(2:3,1), C(2:3,2), C(2:3,3),'o-','LineWidth',5, 'Color', 'b');
plot3(D(2:3,1), D(2:3,2), D(2:3,3),'o-','LineWidth',5, 'Color', 'y');
grid on;

xlabel('X')
ylabel('Y')
zlabel('Z')

% cal shadow length from bottom of object to edges of shadow
s_a = shadowVectorLength(A);
s_b = shadowVectorLength(B);
s_c = shadowVectorLength(C);
s_d = shadowVectorLength(D);

%computing unit vector for ray A, B, C, D
unitVectRay_A = calUnitVector(A(3,:), A(2,:));
unitVectRay_B = calUnitVector(B(3,:), B(2,:));
unitVectRay_C = calUnitVector(C(3,:), C(2,:));
unitVectRay_D = calUnitVector(D(3,:), D(2,:));

% manual tracing by add scalar value
s = 1200;
p1 = A(2,:) + s*unitVectRay_A;
p2 = B(2,:) + s*unitVectRay_B;
p3 = C(2,:) + s*unitVectRay_C;
p4 = D(2,:) + s*unitVectRay_D;

% plot line with manual tracing
plot3([A(2,1) p1(1)], [A(2,2) p1(2)], [A(2,3) p1(3)], '-', 'LineWidth', 1);
plot3([B(2,1) p2(1)], [B(2,2) p2(2)], [B(2,3) p2(3)], '-', 'LineWidth', 1);
plot3([C(2,1) p3(1)], [C(2,2) p3(2)], [C(2,3) p3(3)], '-', 'LineWidth', 1);
plot3([D(2,1) p4(1)], [D(2,2) p4(2)], [D(2,3) p4(3)], '-', 'LineWidth', 1);


% compute intersect position via 2-rays tracing
p_IT_ab = rayintersect(A(2,:), B(2,:), unitVectRay_A, unitVectRay_B);
p_IT_bc = rayintersect(B(2,:), C(2,:), unitVectRay_B, unitVectRay_C);
p_IT_cd = rayintersect(C(2,:), D(2,:), unitVectRay_C, unitVectRay_D);
p_IT_ad = rayintersect(A(2,:), D(2,:), unitVectRay_A, unitVectRay_D);

% ploting 4 intersection points
plot3(p_IT_ab(1), p_IT_ab(2), p_IT_ab(3), '*', 'LineWidth', 5, 'Color', 'r');
plot3(p_IT_bc(1), p_IT_bc(2), p_IT_bc(3), '*', 'LineWidth', 5, 'Color', 'g');
plot3(p_IT_cd(1), p_IT_cd(2), p_IT_cd(3), '*', 'LineWidth', 5, 'Color', 'b');
plot3(p_IT_ad(1), p_IT_ad(2), p_IT_ad(3), '*', 'LineWidth', 5, 'Color', 'y');


% Vertical ray angle figure plot
figure('Name','Elevation ray');
hold on;
view(0,0)  % YZ View
plot3([B(2,1) p_IT_bc(1) C(2,1)], [B(2,2) p_IT_bc(2) C(2,2)], [B(2,3) p_IT_bc(3) C(2,3)], '*-', 'LineWidth', 1, 'Color', 'r');
plot3([A(2,1) p_IT_ad(1) D(2,1)], [A(2,2) p_IT_ad(2) D(2,2)], [A(2,3) p_IT_ad(3) D(2,3)], '*-', 'LineWidth', 1, 'Color', 'g');
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')


% Horizontal ray angle figure plot
figure('Name','Azimuth ray');
view(90,0)   % XZ View
hold on;
plot3([A(2,1) p_IT_ab(1) B(2,1)], [A(2,2) p_IT_ab(2) B(2,2)], [A(2,3) p_IT_ab(3) B(2,3)], '*-', 'LineWidth', 1, 'Color', 'b');
plot3([C(2,1) p_IT_cd(1) D(2,1)], [C(2,2) p_IT_cd(2) D(2,2)], [C(2,3) p_IT_cd(3) D(2,3)], '*-', 'LineWidth', 1, 'Color', 'm');
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')

virLightPos = (p_IT_ab + p_IT_cd)/2;

% azimuth triangle angle calculation
az_u = computeMag(virLightPos(1:2), C(2,1:2));
az_v = computeMag(virLightPos(1:2), D(2,1:2));
az_w = computeMag(C(2,1:2), D(2,1:2));
azimuth_a = acos((az_u^2 + az_w^2 - az_v^2)/(2 * az_u * az_w));
azimuth_b = acos((az_w^2 + az_v^2 - az_u^2)/(2 * az_w * az_v));
azimuth_c = acos((az_v^2 + az_u^2 - az_w^2)/(2 * az_v * az_u));


%2d plot (top view) of azimuth angle
figure('Name','2D top view azimuth angle');
plot([C(2,1) virLightPos(1) D(2,1)], [D(2,2) virLightPos(2) C(2,2)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')

% elevation triangle angle calculation (*Note: in real world unit use x coordinate than y in image coordinate)
ev_u = computeMag([virLightPos(1) virLightPos(3)], [D(2,1) D(2,3)]);
ev_v = computeMag([virLightPos(1) virLightPos(3)], [A(2,1) A(2,3)]);
ev_w = computeMag([D(2,1) D(2,1)], [A(2,3) A(2,3)]);
elevation_a = acos((ev_u^2 + ev_w^2 - ev_v^2)/(2 * ev_u * ev_w));
elevation_b = acos((ev_w^2 + ev_v^2 - ev_u^2)/(2 * ev_w * ev_v));
elevation_c = acos((ev_v^2 + ev_u^2 - ev_w^2)/(2 * ev_v * ev_u));

%2d plot (side view) of elevation angle
figure('Name','2D side view elevation angle');
plot([A(2,1) virLightPos(1) D(2,1)], [A(2,3) virLightPos(3) D(2,3)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Z')

    % nearing the x coordinat make them perfectly square
    function exact = nearestPointY(origin,compare)
        exact = compare;
        diff = abs(origin - compare);
        if(diff(2,1) < 10) 
            exact(2,1) = origin(2,1);
        else
            exact(2,1) = compare(2,1);
        end
            
        if(diff(2,2) < 10) 
            exact(2,2) = origin(2,2);
        else
            exact(2,2) = compare(2,2);
        end
    end
    
    % nearing the y coordinate make them perfectly square
    function exact = nearestPointX(origin,compare)
        exact = compare;
        diff = abs(origin - compare);
        if(diff(1,1) < 10) 
            exact(1,1) = origin(1,1);
        else
            exact(1,1) = compare(1,1);
        end
            
        if(diff(1,2) < 10) 
            exact(1,2) = origin(1,2);
        else
            exact(1,2) = compare(1,2);
        end
    end
    
    % transform to world coordinate
    function worldUnit = homographyTransform(homographyMatrix, inputMatrix)
        PR = homographyMatrix * inputMatrix;
        newX = PR(1,1)/PR(3,1);
        newY = PR(2,1)/PR(3,1);
        worldUnit = [newX; newY];
    end
    
    % function for calculate shadow length
    function ret = shadowVectorLength(position)
        ret = [abs(position(1,1)-position(3,1)), abs(position(1,2)-position(3,2))];
    end

    
    %function for unit vector calculation 
    function unitVect = calUnitVector(origin, dest)
        %l1 = [abs(origin(1) - dest(1)); abs(origin(2) - dest(2)); abs(origin(3) - dest(3))];
        l1 = [origin(1) - dest(1); origin(2) - dest(2); origin(3) - dest(3)];
        mag = sqrt((l1(1) * l1(1)) + (l1(2) * l1(2)) + (l1(3) * l1(3)));
    %     unitVect = [origin(1) + (l1(1)/mag); origin(2) + (l1(2)/mag)];
        unitVect = [l1(1)/mag, l1(2)/mag, l1(3)/mag];
    end

    % function for vector magnitude calculation
    function magValue = computeMag(origin, dest)
        diffVect= origin - dest;
        magValue = sqrt((diffVect(1) * diffVect(1)) + (diffVect(2) * diffVect(2)));
    end

    