function [Height, sliceNumber] = calHeightFromShadow(homographyMatrix, centroidLine, binSegmentShadowImage, virLightPosIMG, virLightPos, interpolateScalar)
n=1;
k=1;
    for i=1:size(centroidLine,1)
        % calculate ray from centroid line to get direction to the lower shadow Edge
        % 1st calculate direction vector (unit vector)
        centroidLine2LightUVect = unitVector2D(double(centroidLine(i,1:2)), virLightPosIMG(1:2)); 
        flag = 0;
        flagFoundSomeEdge = 0;
        shadowLength(n,1) = centroidLine(i,1);
        shadowLength(n,2) = 0;
        Height(n,1) = 0;
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
               % Convert point of shadow edges that respect with virtual light direction
           PtC = [homographyTransform(homographyMatrix, [double(centroidLine(i,1:2)) 1]')' 0];
           Pt0 = [homographyTransform(homographyMatrix, [Pt0 1]')' 0];
           Pt1 = [homographyTransform(homographyMatrix, [Pt1 1]')' 0];
           diffPt = Pt0 - Pt1;
           diffPtC = PtC - Pt1;
           shadowLength = sqrt((diffPt(1) * diffPt(1)) + (diffPt(2) * diffPt(2)));
           shadowPCLength = sqrt((diffPtC(1) * diffPtC(1)) + (diffPtC(2) * diffPtC(2)));

           % find base length between shadow edge and virtual light source (given w)
           uLength = calMag(Pt1, virLightPos);
           vLength = calMag([virLightPos(1:2) 0], virLightPos);
           wLength = calMag(Pt1, [virLightPos(1:2) 0]);   
           angle = acos((uLength^2 + wLength^2 - vLength^2)/(2 * uLength * wLength));
           Height(n,1) = shadowLength * angle;
           HeightPtc(n,1) = shadowPCLength * angle;
           realHeight(n,1) = HeightPtc(n,1) - Height(n,1);
           if(k>50)    
            plot3([virLightPos(1) Pt1(1)], [virLightPos(2) Pt1(2)], [virLightPos(3) Pt1(3)], '*-', 'LineWidth', 1); 
            plot3([virLightPos(1) virLightPos(1)], [virLightPos(2) virLightPos(2)], [virLightPos(3) 0], '*-', 'LineWidth', 1); 
            plot3([virLightPos(1) Pt1(1)], [virLightPos(2) Pt1(2)], [0 0], '*-', 'LineWidth', 1);
            plot3([Pt0(1) Pt1(1)], [Pt0(2) Pt1(2)], [0 0], '*-', 'LineWidth', 1); 
            plot3([Pt0(1) Pt0(1)],[Pt0(2) Pt0(2)],[0 Height(n,1)], 'o-', 'LineWidth', 1);
            k=0;
           end
        end
        n = n + 1;
        k = k + 1;
    end
    sliceNumber = n - 1;
end

