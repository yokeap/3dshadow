clear;
clc;
o1 = [1;5];
o2 = [4;2];
p = [8;10];
% l1 = [abs(o1(1) - p(1)); abs(o1(2) - p(2))];
% mag = sqrt((l1(1) * l1(1)) + (l1(2) * l1(2)));
% unit_o1 = [o1(1) + (l1(1)/mag); o1(2) + (l1(2)/mag)];
unitVect1 = calUnitVector(p, o1);
unitVect2 = calUnitVector(p,o2);
plot([o1(1),p(1)], [o1(2), p(2)], 'Color', 'r', 'LineWidth', 1);
hold on;
plot([o2(1),p(1)], [o2(2), p(2)], 'Color', 'b', 'LineWidth', 1);

uv1 = o1 + unitVect1;
uv2 = o2 + unitVect2;

plot([o1(1),uv1(1)], [o1(2), uv1(2)], 'Color', 'g', 'LineWidth', 3);
plot([o2(1),uv2(1)], [o2(2), uv2(2)], 'Color', 'y', 'LineWidth', 3);

% plot([o1(1),unitVect1(1)], [o1(2), unitVect1(2)], 'Color', 'g', 'LineWidth', 3);
% plot([o2(1),unitVect2(1)], [o2(2), unitVect2(2)], 'Color', 'y', 'LineWidth', 3);

u1 = calScalarLine(o1(1), p(1), unitVect1(1));
v1 = calScalarLine(o1(2), p(2), unitVect1(2));

u2 = calScalarLine(o2(1), p(1), unitVect2(1));
v2 = calScalarLine(o2(2), p(2), unitVect2(2));

%function for unit vector calculation 
function unitVect = calUnitVector(origin, dest)
    %l1 = [abs(origin(1) - dest(1)); abs(origin(2) - dest(2))];
    l1 = [origin(1) - dest(1); origin(2) - dest(2)];
    mag = sqrt((l1(1) * l1(1)) + (l1(2) * l1(2)));
%     unitVect = [origin(1) + (l1(1)/mag); origin(2) + (l1(2)/mag)];
    unitVect = [l1(1)/mag; l1(2)/mag];
end

function scaline = calScalarLine(origin, dest, unitVector)
    scaline = (dest - origin)/unitVector;
end


