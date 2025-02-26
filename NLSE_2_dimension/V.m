function [ V_xy ] = V( x,y)
global omega lambda
if length(x)==1
    V_xy =  omega*(x^2+y^2) + lambda;

else
    V_xy = omega*(x.^2+y.^2)+lambda;
 
end

