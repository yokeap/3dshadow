
P1=homographyTransform(homographyMatrix,A(:,1));
S1=homographyTransform(homographyMatrix,A(:,2));

P2=homographyTransform(homographyMatrix,B(:,1));
S2=homographyTransform(homographyMatrix,B(:,2));

P3=homographyTransform(homographyMatrix,C(:,1));
S3=homographyTransform(homographyMatrix,C(:,2));

P4=homographyTransform(homographyMatrix,D(:,1));
S4=homographyTransform(homographyMatrix,D(:,2));

figure;
plot([P1(1,1) S1(1,1)],[P1(2,1) S1(2,1)],'*-', 'Color', 'r');
hold on;
plot([P2(1,1) S2(1,1)],[P2(2,1) S2(2,1)],'*-', 'Color', 'g');
plot([P3(1,1) S3(1,1)],[P3(2,1) S3(2,1)],'*-', 'Color', 'b');
plot([P4(1,1) S4(1,1)],[P4(2,1) S4(2,1)],'*-', 'Color', 'y');
hold off;


    function worldUnit = homographyTransform(homographyMatrix, inputMatrix)
        PR = homographyMatrix * inputMatrix;
        newX = PR(1,1)/PR(3,1);
        newY = PR(2,1)/PR(3,1);
        worldUnit = [newX; newY];
    end
    

    