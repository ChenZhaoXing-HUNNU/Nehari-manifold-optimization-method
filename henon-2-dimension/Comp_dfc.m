function [u_1,u_2] = Comp_dfc(u)
%%% the Fourier expansion in - \theta
%%% u = \sum_{i=0}^M (u_1m*cos(m\theta) + u_2m*sin(m\theta));
global M theta
len = length(u(1,:));
for m=0:M
     if m== 0
        u_1m = 1/(len)*sum(u,2);
        u_2m = zeros(size(u(:,1)));
        
    else
        u_1m = 2/(len)*(u*cos(m*theta));
        u_2m = 2/(len)*(u*sin(m*theta));
      
    end
    u_1(:,m+1) = u_1m;
    u_2(:,m+1) = u_2m;
 end