function [nuDeg]=prandtlMeyer(Ma,gamma)

a = (gamma+1)/(gamma-1);

firstTerm = sqrt(a)*(atan(sqrt(((Ma^2)-1)/a)));

secondTerm = atan(sqrt((Ma^2)-1));

nu = firstTerm - secondTerm;

nuDeg = nu*180/pi;

end
