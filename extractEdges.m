% extract edges of sample with maximum amplitude and lowest amplitude
% return is (x, ymin, ymax) in image coordinate
% this algorithm is running from 0,0 and scan start in column and row
% innner loop

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
    
    



