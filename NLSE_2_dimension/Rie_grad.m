function RgradE = Rie_grad(u)
%%% compute the Riemannian gradient by solving -\Delta \psi + V \psi= u^3 
psi = Sol_Poisson(u); % solving -\Delta \psi + V \psi= u^3 
%% compute the Riemannian gradient
Hgrad_E = u-psi;
Hgrad_G = 2*u - 4*psi;
alpha = inp(Hgrad_E,Hgrad_G)/inp(Hgrad_G,Hgrad_G) ;
RgradE = Hgrad_E -alpha*Hgrad_G;

