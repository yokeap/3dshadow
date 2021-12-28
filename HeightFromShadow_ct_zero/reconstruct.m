function [EdgesUpperWorld, EdgesLowerWorld, skeletonHeight] = reconstruct(homographyMatrix, binSegmentImage, binSegmentShadowImage, virLightPosIMG, virLightPos, interpolateScalar)
%RECONSTRUCT Summary of this function goes here
%   Detailed explanation goes here

    loop = 1;
    EdgesUpperWorld = [];
    EdgesLowerWorld = [];
    skeletonHeight = [];
    Height = [];
    
    k = 1;
    n = 1;

    % Sample skeletonization
    skBinSegmentImage = bwskel(binSegmentImage);
    figure('Name','Skeleton');
    imshow(skBinSegmentImage)
    
    figure('Name','3D Points Cloud');
    hold on;
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    % 1st scan sample image
    for x = 1:size(binSegmentImage,2)          % scan for x(row)
        flag = 0;
        for y = 1:size(binSegmentImage,1)             % scan for y(column)
            pixValue = binSegmentImage(y,x);
            % edge upper found 
            if flag == 0 && pixValue > 0    
                flag = 1;                   %if found upper set flag
                EdgesUpper(loop,1) = x;
                EdgesUpper(loop,2) = y;
                EdgesUpper(loop,3) = 1;    % add homogeneous coordinate 
                EdgesUpperWorld(loop,:) = homographyTransform(homographyMatrix, [x;y;1])';
                
%                 EdgesLower(loop,1) = x;
%                 EdgesLower(loop,2) = y;
%                 EdgesLower(loop,3) = 1;    % add homogeneous coordinate 
%                 EdgesLowerWorld(loop,:) = homographyTransform(homographyMatrix, [x;y;1])';
                plot3(EdgesUpperWorld(loop,1), EdgesUpperWorld(loop,2), 0 , '.-', 'LineWidth', 1);
            end
            % edge lower found 
            if flag == 1 && pixValue < 1 
                EdgesLower(loop,1) = x;
                EdgesLower(loop,2) = y;
                EdgesLower(loop,3) = 1;    % add homogeneous coordinate 
                EdgesLowerWorld(loop,:) = homographyTransform(homographyMatrix, [x;y;1])';
                plot3(EdgesLowerWorld(loop,1), EdgesLowerWorld(loop,2), 0 , '.-', 'LineWidth', 1);
                Height(loop,1) = 0;
                skeletonHeight(loop,:) = [EdgesLowerWorld(loop,:) 0];
                disp(y)
                disp(y:-1:1)
                for yy = y:-1:1      % scan for y(column)
                    flagHeight = 0;                            % flag for check loop (set when found shadow edge)
                    flagFoundSomeEdge = 0;
                    pixValue = skBinSegmentImage(yy,x);
                    if pixValue > 0
                        skeletonHeight(loop,:) = [homographyTransform(homographyMatrix, [double(x); double(yy); 1])' 0];
                        plot3(skeletonHeight(loop,1), skeletonHeight(loop,2), 0 , '.-', 'LineWidth', 1);
                        % cal unit vector from skeleton position which have
                        % direction relate with light source position
                        skeleton2LightUVect = unitVector2D([double(x), double(yy)], virLightPosIMG(1:2));
                        % interpolate the unit vector.
                        for s=1:interpolateScalar
                            Pt = [double(x), double(yy)] + s*(skeleton2LightUVect);
                            % check value the interpolated position 
                            if flagHeight == 0 && (binSegmentShadowImage(int32(Pt(2)), int32(Pt(1))) > 0)
                                flagHeight = 1;
                                Pt0 = Pt;
                            end
                            if flagHeight == 1 && (binSegmentShadowImage(int32(Pt(2)), int32(Pt(1))) < 1) 
                                flagHeight = 0;
                                Pt1 = Pt; 
                                flagFoundSomeEdge = 1;
                            end
                            if(flagFoundSomeEdge)
                                flagFoundSomeEdge = 0;
                                ptC = [homographyTransform(homographyMatrix, [double(x); double(yy); 1])' 0];
                                ptUpper = [homographyTransform(homographyMatrix, [Pt0 1]')' 0];
                                ptLower = [homographyTransform(homographyMatrix, [Pt1 1]')' 0];
                                diffPt = ptC - ptLower ;
                                shadowLength = sqrt((diffPt(1) * diffPt(1)) + (diffPt(2) * diffPt(2)));
                                % find base length between shadow edge and virtual light source (given w)
                                uLength = calMag(ptLower, virLightPos);
                                vLength = calMag([virLightPos(1:2) 0], virLightPos);
                                wLength = calMag(ptLower, [virLightPos(1:2) 0]);   
                                angle = acos((uLength^2 + wLength^2 - vLength^2)/(2 * uLength * wLength));
                                %Height(loop,1) = shadowLength * angle;
                                skeletonHeight(loop,3) = shadowLength * angle;
                                if(k>50)    
                                    plot3([virLightPos(1) ptLower(1)], [virLightPos(2) ptLower(2)], [virLightPos(3) ptLower(3)], '*-', 'LineWidth', 1); 
                                    plot3([virLightPos(1) virLightPos(1)], [virLightPos(2) virLightPos(2)], [virLightPos(3) 0], '*-', 'LineWidth', 1); 
                                    plot3([virLightPos(1) ptLower(1)], [virLightPos(2) ptLower(2)], [0 0], '*-', 'LineWidth', 1);
                                    plot3([ptUpper(1) ptLower(1)], [ptUpper(2) ptLower(2)], [0 0], '*-', 'LineWidth', 1); 
                                    plot3([ptUpper(1) ptUpper(1)],[ptUpper(2) ptUpper(2)],[0 Height(n,1)], 'o-', 'LineWidth', 1);
                                    k=0;
                                end
                                k = k + 1;
                            end                                                   
                        end              
                    end
                    
                end
                loop = loop + 1;
                flag = 0;                  % if found lower reset flag
            end
            % if found only upper reset flag
            if(y == size(binSegmentImage,1))
                flag = 0;
            end
        end
    end
%     
%     meanVal = mean(skeletonHeight(:, 3));
%     for i=1:size(skeletonHeight(:,3),1)
%         if(skeletonHeight(i,3) < 1)
%             skeletonHeight(i,3) = meanVal/2;
%         end
%     end
%     
%     minVal = min(skeletonHeight(:, 3));
%     for i=1:size(skeletonHeight(:,3),1)
%         if(skeletonHeight(i,3) < 1)
%             skeletonHeight(i,3) = minVal;
%         end
%     end
    
    
% 
%     % finding shadow length from skeleton position and store the shadow
%     % edge at same time.
%     for x = EdgesLower(1,1):EdgesLower(size(EdgesLower,1) ,1) 
%         for y = max(EdgesLower(:,2)):-1:1      % scan for y(column)
%             flag = 0;                            % flag for check loop (set when found shadow edge)
%             flagFoundSomeEdge = 0;
%             pixValue = skBinSegmentImage(y,x);
%             if pixValue > 0
%                 % cal unit vector from skeleton position which have
%                 % direction relate with light source position
%                 skeleton2LightUVect = unitVector2D([double(x), double(y)], virLightPosIMG(1:2));
%                 % interpolate the unit vector.
%                 for s=1:interpolateScalar
%                     Pt = [double(x), double(y)] + s*(skeleton2LightUVect);
%                     % check value the interpolated position 
%                     if flag == 0 && (binSegmentShadowImage(int32(Pt(2)), int32(Pt(1))) > 0)
%                         flag = 1;
%                         Pt0 = Pt;
%                         
%                     end
%                     if flag == 1 && (binSegmentShadowImage(int32(Pt(2)), int32(Pt(1))) < 1) 
%                         flag = 0;
%                         Pt1 = Pt;
%                         flagFoundSomeEdge = 1;
%                         break;
%                     end
%                 end
%                 if(flagFoundSomeEdge)
%                     skeletonSample(n,:) = [homographyTransform(homographyMatrix, [double(x); double(y); 1])' 0];
%                     plot3(skeletonSample(n,1), skeletonSample(n,2), 0 , '.-', 'LineWidth', 1);
%                     ptC = [homographyTransform(homographyMatrix, [double(x); double(y); 1])' 0];
%                     ptUpper = [homographyTransform(homographyMatrix, [Pt0 1]')' 0];
%                     ptLower = [homographyTransform(homographyMatrix, [Pt1 1]')' 0];
% %                     ptUpper = [homographyTransform(homographyMatrix, [Pt0 1]')' 0];
% %                     ptLower = [homographyTransform(homographyMatrix, [Pt1 1]')' 0];
%                     diffPt = ptUpper - ptLower ;
%                     shadowLength = sqrt((diffPt(1) * diffPt(1)) + (diffPt(2) * diffPt(2)));
%                     % find base length between shadow edge and virtual light source (given w)
%                     uLength = calMag(ptLower, virLightPos);
%                     vLength = calMag([virLightPos(1:2) 0], virLightPos);
%                     wLength = calMag(ptLower, [virLightPos(1:2) 0]);   
%                     angle = acos((uLength^2 + wLength^2 - vLength^2)/(2 * uLength * wLength));
%                     Height(n,1) = shadowLength * angle;
%                     if(k>50)    
%                         plot3([virLightPos(1) ptLower(1)], [virLightPos(2) ptLower(2)], [virLightPos(3) ptLower(3)], '*-', 'LineWidth', 1); 
%                         plot3([virLightPos(1) virLightPos(1)], [virLightPos(2) virLightPos(2)], [virLightPos(3) 0], '*-', 'LineWidth', 1); 
%                         plot3([virLightPos(1) ptLower(1)], [virLightPos(2) ptLower(2)], [0 0], '*-', 'LineWidth', 1);
%                         plot3([ptUpper(1) ptLower(1)], [ptUpper(2) ptLower(2)], [0 0], '*-', 'LineWidth', 1); 
%                         plot3([ptUpper(1) ptUpper(1)],[ptUpper(2) ptUpper(2)],[0 Height(n,1)], 'o-', 'LineWidth', 1);
%                         k=0;
%                     end
%                     n = n + 1;
%                     k = k + 1;
%                 end              
%             end
%         end
%     end

end
