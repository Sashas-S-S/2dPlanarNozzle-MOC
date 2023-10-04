function [Ma]=inversePrandtl(nu,gamma)

%nuDeg = nu/180;

f = @(Ma) nu - prandtlMeyer(Ma,gamma)

Ma = fsolve(f,1);

end