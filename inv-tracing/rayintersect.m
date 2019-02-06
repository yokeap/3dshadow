clear;
clc;
o1 = [1;5];
o2 = [4;2];
% unitVect1 = [1.813733471206735,5.581238193719097];
% unitVect2 = [4.447213595499958,2.894427190999916];

unitVect1 = [0.813733471206735;0.581238193719097];
unitVect2 = [0.447213595499958;0.894427190999916];


% %rearrange to Ax = B
% % A = [unitVect1(1) -unitVect2(1); unitVect1(2) -unitVect2(2)];
% % B =[o2(1)-o1(1); o2(2)-o1(2)];
% 
% A = [unitVect1(1) -unitVect2(1)];
% B = [o2(1)-o1(1)];
% X = A/B
% 


W = o1 - o2;

u_dot_u = dot(unitVect1, unitVect1);
v_dot_v = dot(unitVect2, unitVect2);
u_dot_v = dot(unitVect1, unitVect2);
w_dot_u = dot(W, unitVect1);
w_dot_v = dot(W, unitVect2);

denom = (u_dot_u * v_dot_v) - (u_dot_v * u_dot_v);

s =  ((u_dot_v/denom) * w_dot_v) - ((v_dot_v/denom) * w_dot_u);
t = -((u_dot_v/denom) * w_dot_u) + ((u_dot_u/denom) * w_dot_v);

p = (o1 + s*unitVect1) + (o2 + t*unitVect2);



