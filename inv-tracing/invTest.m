clear;
clc;

A = [904, 477, 0;
     904, 248, 0;
     1217, 248, 0];
    
B = [2830, 477, 0;
     2830, 248, 0;
     3333, 248, 0];
 
C = [2830, 2277, 0;
     2830, 2522, 0;
     3333, 2522, 0];
 
D = [904, 2277, 0;
     904, 2522, 0;
     1209, 2522, 0];

A = img2cart(A);
B = img2cart(B);
C = img2cart(C);
D = img2cart(D);

plot3(A(:,1), A(:,2), A(:,3),'o-');
hold on
plot3(B(:,1), B(:,2), B(:,3),'o-');
plot3(C(:,1), C(:,2), C(:,3),'o-');
plot3(D(:,1), D(:,2), D(:,3),'o-');
grid on



function newCoordinate = img2cart(coordinate)
    newCoordinate = coordinate;
    newCoordinate(:,2) = 2648 - coordinate(:,2);
end