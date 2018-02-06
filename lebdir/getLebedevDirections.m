function P = getLebedevDirections(N)
% GETLEBEDEVDIRECTIONS    Get unit vectors well spread out over sphere
%
% P = getLebedevDirections(N) returns N times 3 array with a unit vector in
% each row. Unit vectors are Lebedev quadrature points from the function
% getLebedevSphere, where duplicate vectors up to the sign have been
% removed.
%
% Valid number of directions are: 3, 7, 13, 19, 25, 37, 43, 55, 73, 85, 97,
% 115, 133, 151, 175, 217, 295, 385, 487, 601, 727, 865, 1015, 1177, 1351,
% 1537, 1735, 1945, 2167, 2401, 2647, 2905

leb = getLebedevSphere(2*N);
P = [leb.x,leb.y,leb.z];
P = getHalf(P);