    %function for unit vector calculation 
    function unitVect = calUnitVector(origin, dest)
        %l1 = [abs(origin(1) - dest(1)); abs(origin(2) - dest(2)); abs(origin(3) - dest(3))];
        l1 = [origin(1) - dest(1); origin(2) - dest(2); origin(3) - dest(3)];
        mag = sqrt((l1(1) * l1(1)) + (l1(2) * l1(2)) + (l1(3) * l1(3)));
    %     unitVect = [origin(1) + (l1(1)/mag); origin(2) + (l1(2)/mag)];
        unitVect = [l1(1)/mag, l1(2)/mag, l1(3)/mag];
    end
