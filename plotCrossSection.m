function out = plotCrossSection(ptClound, section)

    %New figure establishing
    figure
    %100th slicing plot
    plot(ptClound(section,4:5), ptClound(section,2), '*')
    hold on
    plot(ptClound(section,3), ptClound(section,6), '*') 
    plot(ptClound(section,3), ptClound(section,7), '*')
    %center of slice
    plot(ptClound(section,3), ptClound(section,2), '*')

end

