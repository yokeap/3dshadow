% function for vector magnitude calculation
function magValue = computeMag(origin, dest)
    diffVect= origin - dest;
    magValue = sqrt((diffVect(1) * diffVect(1)) + (diffVect(2) * diffVect(2)));
end