function NewVect = calRotateYAxis(uuVect, theta)
%CALROTATEXAXIS Summary of this function goes here
%   Detailed explanation goes here
%     NewVect =[1 0 0;0 cosd(theta) -sind(theta);0 sind(theta) cosd(theta)] * uuVect';
    NewVect =[cosd(theta) 0 sin(theta);0  1 0;-sind(theta) 0 cosd(theta)] * uuVect';
end

