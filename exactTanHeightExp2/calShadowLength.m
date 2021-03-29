function [heightWorldUnit, sliceNumber] = calShadowLength(centroidLine, binSegmentShadowImage, virLightPos, interpolateScalar)
n=1;
    for i=1:size(centroidLine,1)
        % calculate ray from centroid line to get direction to the lower shadow Edge
        % 1st calculate direction vector (unit vector)
        centroidLine2LightUVect = unitVector2D(double(centroidLine(i,1:2)), virLightPos(1:2)); 
        flag = 0;
        flagFoundSomeEdge = 0;
        HeightWorldUnit(n,1) = centroidLine(i,1);
        HeightWorldUnit(n,2) = 0;
        for s=1:interpolateScalar
            Pt = double(centroidLine(i,1:2)) + s*(centroidLine2LightUVect);
            if flag == 0 && (binSegmentShadowImage(int32(Pt(2)), int32(Pt(1))) > 0)
                flag = 1;
                Pt0 = Pt;
                flagFoundSomeEdge = 1;
            end
            if flag == 1 && (binSegmentShadowImage(int32(Pt(2)), int32(Pt(1))) < 1) 
                flag = 0;
                Pt1 = Pt;
                break;
            end
        end
        if(flagFoundSomeEdge)
           HeightWorldUnit(n,2) = sqrt((Pt0(1) * Pt0(1)) + (Pt1(2) * Pt1(2)));
            plot3([Pt0(1) Pt1(1)], [Pt0(2) Pt1(2)], [0 0], '*-', 'LineWidth', 1); 
        end
        n = n + 1;
    end
    sliceNumber = n - 1;
end

