clear;
clc;
close all;

img_max = [3514, 2648];

A = [904, 477, 0;
     1215, 399, 0;
     904, 477, 800];
    
B = [2850, 477, 0;
     3336, 399, 0;
     2850, 477, 800];
 
C = [2850, 2277, 0;
     3336, 2377, 0;
     2850, 2277, 800];
 
D = [904, 2277, 0;
     1215, 2377, 0;
     904, 2277, 800];

A = img2cart(A, img_max);
B = img2cart(B, img_max);
C = img2cart(C, img_max);
D = img2cart(D, img_max);

% plot3(A(:,1), A(:,2), A(:,3),'o-','LineWidth',2);
% hold on
% plot3(B(:,1), B(:,2), B(:,3),'o-','LineWidth',2);
% plot3(C(:,1), C(:,2), C(:,3),'o-','LineWidth',2);
% plot3(D(:,1), D(:,2), D(:,3),'o-','LineWidth',2);
% grid on

% plot vector from shadow edges to top of object
figure('Name','Inverese Ray Tracing');
hold on;
plot3(A(2:3,1), A(2:3,2), A(2:3,3), 'o-','LineWidth',5);
plot3(B(2:3,1), B(2:3,2), B(2:3,3),'o-','LineWidth',5);
plot3(C(2:3,1), C(2:3,2), C(2:3,3),'o-','LineWidth',5);
plot3(D(2:3,1), D(2:3,2), D(2:3,3),'o-','LineWidth',5);
grid on;

xlabel('X')
ylabel('Y')
zlabel('Z')

% cal shadow length from bottom of object to edges of shadow
s_a = shadowVectorLength(A);
s_b = shadowVectorLength(B);
s_c = shadowVectorLength(C);
s_d = shadowVectorLength(D);

%plot3(A(1,1), A(3,1), 0,'o-','LineWidth',2);

%computing unit vector for ray A, B, C, D
unitVectRay_A = calUnitVector(A(3,:), A(2,:));
unitVectRay_B = calUnitVector(B(3,:), B(2,:));
unitVectRay_C = calUnitVector(C(3,:), C(2,:));
unitVectRay_D = calUnitVector(D(3,:), D(2,:));

% plot3(unitVectRay_A(1), unitVectRay_A(2), unitVectRay_A(3),'o-','LineWidth',2);
% hold on
% plot3(unitVectRay_B(1), unitVectRay_B(2), unitVectRay_B(3),'o-','LineWidth',2);
% plot3(unitVectRay_C(1), unitVectRay_C(2), unitVectRay_C(3),'o-','LineWidth',2);
% plot3(unitVectRay_D(1), unitVectRay_D(2), unitVectRay_D(3),'o-','LineWidth',2);
% grid on

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
xlabel('X')
ylabel('Y')
zlabel('Z')

% compute intersect position via 2-rays tracing
p_IT_ab = rayintersect(A(2,:), B(2,:), unitVectRay_A, unitVectRay_B);
p_IT_bc = rayintersect(B(2,:), C(2,:), unitVectRay_B, unitVectRay_C);
p_IT_cd = rayintersect(C(2,:), D(2,:), unitVectRay_C, unitVectRay_D);
p_IT_ad = rayintersect(A(2,:), D(2,:), unitVectRay_A, unitVectRay_D);
xlabel('X')
ylabel('Y')
zlabel('Z')

% p = A(2,:) + unitVectRay_A;

% ploting 4 intersection points
plot3(p_IT_ab(1), p_IT_ab(2), p_IT_ab(3), '*', 'LineWidth', 5);
plot3(p_IT_bc(1), p_IT_bc(2), p_IT_bc(3), '*', 'LineWidth', 5);
plot3(p_IT_cd(1), p_IT_cd(2), p_IT_cd(3), '*', 'LineWidth', 5);
plot3(p_IT_ad(1), p_IT_ad(2), p_IT_ad(3), '*', 'LineWidth', 5);
xlabel('X')
ylabel('Y')
zlabel('Z')


% Horizontal ray angle figure plot
figure('Name','Azimuth ray');
hold on;
view(90,0)  % YZ View
plot3([B(2,1) p_IT_bc(1) C(2,1)], [B(2,2) p_IT_bc(2) C(2,2)], [B(2,3) p_IT_bc(3) C(2,3)], '*-', 'LineWidth', 1);
plot3([A(2,1) p_IT_ad(1) D(2,1)], [A(2,2) p_IT_ad(2) D(2,2)], [A(2,3) p_IT_ad(3) D(2,3)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')


% vertical ray angle figure plot
figure('Name','Elevation ray');
view(0,0)   % XZ View
hold on;
plot3([A(2,1) p_IT_ab(1) B(2,1)], [A(2,2) p_IT_ab(2) B(2,2)], [A(2,3) p_IT_ab(3) B(2,3)], '*-', 'LineWidth', 1);
plot3([C(2,1) p_IT_cd(1) D(2,1)], [C(2,2) p_IT_cd(2) D(2,2)], [C(2,3) p_IT_cd(3) D(2,3)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')

% new figure plot, just 2 sides: bc, cd for preparing finding the base
% anglel
figure('Name','Virtual Light source radiation azimuth and altitude angles');
hold on
plot3([B(2,1) p_IT_bc(1) C(2,1)], [B(2,2) p_IT_bc(2) C(2,2)], [B(2,3) p_IT_bc(3) C(2,3)], '*-', 'LineWidth', 1);
plot3([A(2,1) p_IT_ad(1) D(2,1)], [A(2,2) p_IT_ad(2) D(2,2)], [A(2,3) p_IT_ad(3) D(2,3)], '*-', 'LineWidth', 1);
plot3([C(2,1) p_IT_cd(1) D(2,1)], [C(2,2) p_IT_cd(2) D(2,2)], [C(2,3) p_IT_cd(3) D(2,3)], '*-', 'LineWidth', 1);
plot3([A(2,1) p_IT_ab(1) B(2,1)], [A(2,2) p_IT_ab(2) B(2,2)], [A(2,3) p_IT_ab(3) B(2,3)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')

% azimuth triangle angle calculation
az_u = computeMag(p_IT_bc(1,1:2), B(2,1:2));
az_v = computeMag(B(2,1:2), C(2,1:2));
az_w = computeMag(C(2,1:2), p_IT_bc(1,1:2));
azimuth_a = acosd((az_u^2 + az_w^2 - az_v^2)/(2 * az_u * az_w));
azimuth_b = acosd((az_w^2 + az_v^2 - az_u^2)/(2 * az_w * az_v));
azimuth_c = acosd((az_v^2 + az_u^2 - az_w^2)/(2 * az_v * az_u));

%2d plot (top view) of azimuth angle
figure('Name','2D top view azimuth angle');
plot([B(2,1) p_IT_bc(1) C(2,1)], [B(2,2) p_IT_bc(2) C(2,2)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')

% elevation triangle angle calculation
ev_u = computeMag(p_IT_cd(1,1:2), C(2,1:2));
ev_v = computeMag(C(2,1:2), D(2,1:2));
ev_w = computeMag(D(2,1:2), p_IT_cd(1,1:2));
elevation_a = acosd((ev_u^2 + ev_w^2 - ev_v^2)/(2 * ev_u * ev_w));
elevation_b = acosd((ev_w^2 + ev_v^2 - ev_u^2)/(2 * ev_w * ev_v));
elevation_c = acosd((ev_v^2 + ev_u^2 - ev_w^2)/(2 * ev_v * ev_u));

%2d plot (side view) of elevation angle
figure('Name','2D side view elevation angle');
plot([C(2,1) p_IT_cd(1) D(2,1)], [C(2,2) p_IT_cd(2) D(2,2)], '*-', 'LineWidth', 1);
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')




% function for convert image unit to cartesian coordinate
function newCoordinate = img2cart(coordinate, img_max)
    newCoordinate = coordinate;
    newCoordinate(:,2) = img_max(2) - coordinate(:,2);
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