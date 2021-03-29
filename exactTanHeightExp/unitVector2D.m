function unitVect = unitVector2D(origin, dest)
        l1 = [origin(1) - dest(1); origin(2) - dest(2)];
        mag = sqrt((l1(1) * l1(1)) + (l1(2) * l1(2)));
        unitVect = [l1(1)/mag, l1(2)/mag];
end

