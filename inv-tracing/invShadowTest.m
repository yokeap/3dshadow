clear;
clc;

A = [904, 477, 0;
     1215, 399, 0;
     904, 477, 800];
    
B = [2850, 472, 0;
     3336, 399, 0;
     2850, 472, 800];
 
C = [2850, 2277, 0;
     3336, 2377, 0;
     2850, 2277, 800];
 
D = [904, 2269, 0;
     1215, 2377, 0;
     904, 2269, 800];

A = img2cart(A);
B = img2cart(B);
C = img2cart(C);
D = img2cart(D);

% plot3(A(:,1), A(:,2), A(:,3),'o-','LineWidth',2);
% hold on
% plot3(B(:,1), B(:,2), B(:,3),'o-','LineWidth',2);
% plot3(C(:,1), C(:,2), C(:,3),'o-','LineWidth',2);
% plot3(D(:,1), D(:,2), D(:,3),'o-','LineWidth',2);
% grid on

plot3(A(2:3,1), A(2:3,2), A(2:3,3), 'o-','LineWidth',5);
hold on
plot3(B(2:3,1), B(2:3,2), B(2:3,3),'o-','LineWidth',5);
plot3(C(2:3,1), C(2:3,2), C(2:3,3),'o-','LineWidth',5);
plot3(D(2:3,1), D(2:3,2), D(2:3,3),'o-','LineWidth',5);
grid on

% xlabel('X')
% ylabel('Y')
% zlabel('Z')

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

s = 11000;
p1 = A(2,:) + s*unitVectRay_A;
p2 = B(2,:) + s*unitVectRay_B;
p3 = C(2,:) + s*unitVectRay_C;
p4 = D(2,:) + s*unitVectRay_D;

% p = A(2,:) + unitVectRay_A;
 
plot3([A(2,1) p1(1)], [A(2,2) p1(2)], [A(2,3) p1(3)], '-', 'LineWidth', 2);
plot3([B(2,1) p2(1)], [B(2,2) p2(2)], [B(2,3) p2(3)], '-', 'LineWidth', 2);
plot3([C(2,1) p3(1)], [C(2,2) p3(2)], [C(2,3) p3(3)], '-', 'LineWidth', 2);
plot3([D(2,1) p4(1)], [D(2,2) p4(2)], [D(2,3) p4(3)], '-', 'LineWidth', 2);

function newCoordinate = img2cart(coordinate)
    newCoordinate = coordinate;
    newCoordinate(:,2) = 2648 - coordinate(:,2);
end

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
