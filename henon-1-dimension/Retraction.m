function [result_u,laplace_result_u,result_energy] = Retraction(alpha_0, RgradE,u_k, laplace_u_k, laplace_RgradE)
%%% retraction map: u_{k+1}= R_{u_k}( - alpha_0*RgradE) 
%%% laplace_u_k: \Delta u_k
%%% RgradE: \nabla_{\mathcal{N}} E(u_k)
%%%  laplace_RgradE £º \Delta(\nabla_{\mathcal{N}} E(u_k))
%%% alpha_0 : step-size
global g_xy p h
u_new = u_k - alpha_0*RgradE;
laplace_u_new = laplace_u_k - alpha_0*laplace_RgradE;
up =  h*(-laplace_u_new'*u_new);
down =h*sum(g_xy.*abs(u_new).^(p+1));
rho = (up/down)^(1/(p-1));
result_u = rho*u_new;
laplace_result_u = rho*laplace_u_new;
result_energy = h*(0.5*(-laplace_result_u'*result_u) - 1/(p+1)*sum(g_xy.*abs(result_u).^(p+1)));
end

