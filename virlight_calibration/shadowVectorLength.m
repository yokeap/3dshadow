% function for calculate shadow length
function ret = shadowVectorLength(position)
    ret = [abs(position(1,1)-position(3,1)), abs(position(1,2)-position(3,2))];
end
