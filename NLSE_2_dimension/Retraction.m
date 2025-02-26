function [u_now,energy_now] = Retraction(alpha_0, gradE,u_k)
%%% Nehari retraction map
global h1 h2 
s = h1*h2;
u_new = u_k + alpha_0*gradE;
down = sum((u_new).^(4))*s;
up =inp(u_new,u_new);
sigma = (up/down)^(1/2); 
u_now = sigma*u_new;
energy_now = 0.5*inp(u_now,u_now) - 1/4*sum((u_now).^(4))*s ;
end

