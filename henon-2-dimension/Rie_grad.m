function RgradE = Rie_grad(u)
%%% compute the Riemannian gradient by solving -\Delta \psi = g|u|^{p-1}^u 
global R p l
g_l = R.^l;
f = g_l.*abs(u.^(p-1)).*u;


%%% solve -\Delta \psi(x,y) = f 
psi = Sol_Poisson(f);
%% compute the Riemannian gradient
Hgrad_E = u-psi;
Hgrad_G = 2*u - (p+1)*psi;
alpha = inp(Hgrad_E,Hgrad_G)/inp(Hgrad_G,Hgrad_G) ;
RgradE = Hgrad_E -alpha*Hgrad_G;


