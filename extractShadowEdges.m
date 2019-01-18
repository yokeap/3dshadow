function EdgesCoordinate = extractShadowEdges(binShadowImage, mangoPosition)

loop = 1;

for x = 1:length(mangoPosition)         
    flag = 0;
    for y = 1:size(binShadowImage,1)   
        pixValue = binShadowImage(y,mangoPosition(x));
        if flag == 0 && pixValue > 0
            flag = 1;
            EdgesCoordinate(loop,1) = mangoPosition(x);
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
    
    



