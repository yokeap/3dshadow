
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

input_A = [1.591619136460555e+03,1.426180970149254e+03;7.014912046908303e+02,9.517694562899778e+02];
input_B = [3.873817430703625e+03,4.005319562899787e+03;7.014912046908303e+02,9.517694562899778e+02];
input_C = [3.873817430703625e+03,3.992593550106611e+03;3.174579690831556e+03,3.700588219616204e+03];
input_D = [1.570409115138593e+03,1.413454957356077e+03;3.174579690831556e+03,3.713314232409381e+03];

% homographyMatrix = ...
%        [5.079053647431133e-06,-5.475219832828209e-04,0.224511299503042;
%        -5.494275672572568e-04,-1.011837140387550e-06,0.974457649498132;
%        1.174112215996749e-08,-1.921899063852394e-08,-0.005134593677710];

% % load instrincsic camera calibrated parameters
% load('../calibration/cameraParams.mat');
% % load computed homography matrix
% 
% % load distorted image
% imOrig = imread('../SampleImages/Fish Test/01-03-2019.JPG');
% % get undistorting image and offset origin
% [im, newOrigin] = undistortImage(imOrig, cameraParams, 'OutputView', 'full');
% figure; imshow(im); title('Undistorted Image');
% axis image;
% axis ij ;

% get the coordinate of shadow ray
% hold on;
% input_A = ginput(2)';
% input_B = ginput(2)'; 
% input_C = ginput(2)';
% input_D = ginput(2)';

input_B = nearestPointY(input_A,input_B);
input_C = nearestPointX(input_B,input_C);
input_D = nearestPointX(input_A,input_D);
input_D = nearestPointY(input_C,input_D);

plot(input_A(1,1:2), input_A(2, 1:2), '*-');
hold on;
plot(input_B(1,1:2), input_B(2, 1:2), '*-');
plot(input_C(1,1:2), input_C(2, 1:2), '*-');
plot(input_D(1,1:2), input_D(2, 1:2), '*-');
hold off;

% rearrange the matrix
A = [input_A input_A(:,1); [0 0 52]]';
B = [input_B input_B(:,1); [0 0 52]]';
C = [input_C input_C(:,1); [0 0 52]]';
D = [input_D input_D(:,1); [0 0 52]]';


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
s = 12000;
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
view(90,0)  % YZ View
plot3([B(2,1) p_IT_bc(1) C(2,1)], [B(2,2) p_IT_bc(2) C(2,2)], [B(2,3) p_IT_bc(3) C(2,3)], '*-', 'LineWidth', 1, 'Color', 'r');
plot3([A(2,1) p_IT_ad(1) D(2,1)], [A(2,2) p_IT_ad(2) D(2,2)], [A(2,3) p_IT_ad(3) D(2,3)], '*-', 'LineWidth', 1, 'Color', 'g');
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')


% Horizontal ray angle figure plot
figure('Name','Azimuth ray');
view(0,0)   % XZ View
hold on;
plot3([A(2,1) p_IT_ab(1) B(2,1)], [A(2,2) p_IT_ab(2) B(2,2)], [A(2,3) p_IT_ab(3) B(2,3)], '*-', 'LineWidth', 1, 'Color', 'b');
plot3([C(2,1) p_IT_cd(1) D(2,1)], [C(2,2) p_IT_cd(2) D(2,2)], [C(2,3) p_IT_cd(3) D(2,3)], '*-', 'LineWidth', 1, 'Color', 'm');
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')

virLightPosIMG = (p_IT_cd + p_IT_ab)/2;
virLightPos = [homographyTransform(homographyMatrix, [virLightPosIMG(1); virLightPosIMG(2); 1])' virLightPosIMG(3)];

% azimuth triangle angle calculation
az_u = computeMag(virLightPosIMG(1:2), C(2,1:2));
az_v = computeMag(virLightPosIMG(1:2), D(2,1:2));
az_w = computeMag(C(2,1:2), D(2,1:2));
azimuth_a = acos((az_u^2 + az_w^2 - az_v^2)/(2 * az_u * az_w));
azimuth_b = acos((az_v^2 + az_w^2 - az_u^2)/(2 * az_v * az_w));
azimuth_c = acos((az_u^2 + az_v^2 - az_w^2)/(2 * az_v * az_u));


%2d plot (top view) of azimuth angle
figure('Name','2D top view azimuth angle');
plot([C(2,1) virLightPosIMG(1) D(2,1)], [D(2,2) virLightPosIMG(2) C(2,2)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')

% elevation triangle angle calculation
ev_u = computeMag([virLightPosIMG(2) virLightPosIMG(3)], [D(2,2) D(2,3)]);
ev_v = computeMag([virLightPosIMG(2) virLightPosIMG(3)], [A(2,2) A(2,3)]);
ev_w = computeMag([D(2,2) D(2,3)], [A(2,2) A(2,3)]);
elevation_a = acos((ev_u^2 + ev_w^2 - ev_v^2)/(2 * ev_u * ev_w));
elevation_b = acos((ev_v^2 + ev_w^2 - ev_u^2)/(2 * ev_v * ev_w));
elevation_c = acos((ev_u^2 + ev_v^2 - ev_w^2)/(2 * ev_v * ev_u));

%2d plot (side view) of elevation angle
figure('Name','2D side view elevation angle');
plot([A(2,2) virLightPosIMG(2) D(2,2)], [A(2,3) virLightPosIMG(3) D(2,3)], '*-', 'LineWidth', 1);
grid on;
xlabel('Y')
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

    