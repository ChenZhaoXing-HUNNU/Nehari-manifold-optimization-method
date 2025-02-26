function inp_u_v = inp(u,v)
%%% compute the inner product (u,v) = \int_{\Omega} \nabla u \cdot \nabla dxdy
global  K h1 h2 V_xy
s =h1*h2;
inp_u_v = u'*K*v + sum(V_xy.*u.*v)*s;