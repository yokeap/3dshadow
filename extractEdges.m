function EdgesCoordinate = extractEdges(binShadowImage)

loop = 1;

for x = 1:size(binShadowImage,2)                      
    flag = 0;
    for y = 1:size(binShadowImage,1)   
        pixValue = binShadowImage(y,x);
        if flag == 0 && pixValue > 0
            flag = 1;
            EdgesCoordinate(loop,1) = x;
            EdgesCoordinate(loop,2) = y;
        end
        if flag == 1 && pixValue < 1 
            EdgesCoordinate(loop,3) = y;
            flag = 0;
            loop = loop + 1;
            %fprintf("%i\n", loop);
            break;
        end
    end
end
    
    



