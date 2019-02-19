function EdgesCoordinate = extractEdges(binImage)

loop = 1;

for x = 1:size(binImage,2)                      
    flag = 0;
    for y = 1:size(binImage,1)   
        pixValue = binImage(y,x);
        if flag == 0 && pixValue > 0
            flag = 1;
            EdgesCoordinate(loop,1) = x;
            EdgesCoordinate(loop,2) = y;
            EdgesCoordinate(loop,4) = 1;    % add homogeneous coordinate 
        end
        if flag == 1 && pixValue < 1 
            EdgesCoordinate(loop,3) = y;
            EdgesCoordinate(loop,4) = 1;    % add homogeneous coordinate 
            flag = 0;
            loop = loop + 1;
            %fprintf("%i\n", loop);
            break;
        end
    end
end
    
    



