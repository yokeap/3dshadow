function [position] = rayintersect(origin1, origin2, unitVect1, unitVect2)
    %RAYINTERSECT Summary of this function goes here
    %   compute the position of 2-rays intersection
    %   Author: Siwakorn Sukprasertchai

    W = origin1 - origin2;

    u_dot_u = dot(unitVect1, unitVect1);
    v_dot_v = dot(unitVect2, unitVect2);
    u_dot_v = dot(unitVect1, unitVect2);
    w_dot_u = dot(W, unitVect1);
    w_dot_v = dot(W, unitVect2);

    denom = (u_dot_u * v_dot_v) - (u_dot_v * u_dot_v);

    s =  ((u_dot_v/denom) * w_dot_v) - ((v_dot_v/denom) * w_dot_u);
    t = -((u_dot_v/denom) * w_dot_u) + ((u_dot_u/denom) * w_dot_v);

    position = (origin1 + s*unitVect1) + (origin2 + t*unitVect2);
    
    position = position/2;

end

