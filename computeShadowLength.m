function shadowLength = computeShadowLength(shadowEdges)
for x = 1:size(shadowEdges)
    shadowLength(x,1) = shadowEdges(x);
    shadowLength(x,2) = shadowEdges(x,3) - shadowEdges(x,2);
end
