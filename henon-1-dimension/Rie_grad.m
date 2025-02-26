function [RgradE,laplace_RgradE] = Rie_grad(u, laplace_u)
%%%$ compute the Riemannian gradient \nabla_{N} E(u)
%%% solve poisson方程 -w'' = g|u|^{p-1}u; 
%%% gradE = u-w 
%%% gradG= 2u-(p+1)w
%%% RgradE = gradE - (gradE,gradG)/(gradG,gradG)* gradG...
global g_xy p h K_BD M_cBD

%%% 
 f = g_xy.*abs(u).^(p-1).*u;
F = M_cBD*f; % load vector
omega = K_BD\F; 
% omega = Sol_Poisson(u);
gradE = u - omega;
laplace_gradE = laplace_u + f;
gradG = 2*u - (p+1)*omega;
laplace_gradG = 2*laplace_u + (p+1)*f;
in_EG = -h* laplace_gradE'* gradG; % (gradE,gradG)_H；
in_GG = -h* laplace_gradG'* gradG; % (gradE,gradG)_H之间的内积；
rho  = in_EG/in_GG;
RgradE= gradE - rho*gradG; 
laplace_RgradE = laplace_gradE - rho*laplace_gradG; 
