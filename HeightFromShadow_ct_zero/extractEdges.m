function [EdgesUpper, EdgesLower] = extractEdges(binImage)

loop = 1;

for x = 1:size(binImage,2)                      
    flag = 0;
    for y = 1:size(binImage,1)   
        pixValue = binImage(y,x);
        if flag == 0 && pixValue > 0    
            flag = 1;                   %if found upper set flag
            EdgesUpper(loop,1) = x;
            EdgesUpper(loop,2) = y;
            EdgesUpper(loop,3) = 1;    % add homogeneous coordinate 
        end
        if flag == 1 && pixValue < 1 
            EdgesLower(loop,1) = x;
            EdgesLower(loop,2) = y;
            EdgesLower(loop,3) = 1;    % add homogeneous coordinate 
            flag = 0;                  % if found lower reset flag
            loop = loop + 1;
            %fprintf("%i\n", loop);
            break;
        end
        % if found only upper reset flag
        if(y == size(binImage,1))
            flag = 0;
        end
    end
end
    
    



