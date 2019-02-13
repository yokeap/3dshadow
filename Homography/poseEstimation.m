% R =[0.999754876586451,-0.002339715563935,0.022016186606441;
%     0.003634721363716,0.998253387993276,-0.058965771092672;
%     -0.021839769738207,0.059031339905385,0.998017196929370];
% 
% T = [-1.568660138072978e+02,-67.838402390008700,1.068218666557448e+03];


dat = ginput(4);
datw = [1 1 1 1]';

Data=[dat datw]';
aplot(Data, 'g+')
w=50 ; %mm 5 cm
h=50 ; %mm 5 cm

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

%Measure the coins diameter 1
p=ginput(2)
pw=[1 1]';
P=[p pw]'

PR=H*P
%Computing the distance of diameter of the coins
delx=((PR(1,1)/PR(3,1))-(PR(1,2)/PR(3,2)));
dely=((PR(2,1)/PR(3,1))-(PR(2,2)/PR(3,2)));
Lenght = sqrt(delx^2+dely^2)
aplot(P,'g.-')
text(p(1,1),p(1,2),[num2str(Lenght)]);