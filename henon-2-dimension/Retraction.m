function [ret_res, energy] = Retraction(alpha,eta,u_0)
global  p  g_l

u_new= u_0+alpha*eta;
ret_up = inp(u_new,u_new);
ret_down = integ(g_l.*(abs(u_new).^(p+1)));  
rho_unew = (ret_up/ret_down)^(1/(p-1));
ret_res = rho_unew*u_new;
energy = 0.5*inp(ret_res,ret_res)-1/(p+1)*integ(g_l.*(abs(ret_res).^(p+1)));