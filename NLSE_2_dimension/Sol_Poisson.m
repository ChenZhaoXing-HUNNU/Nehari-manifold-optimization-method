function [psi]=Sol_Poisson(u)
%%% compute the poisson equation -\Delta \psi + V(x,y)\psi = u^3 , (x,y)\in \Omge = (-1,1)^2 
%%% \psi = 0, (x,y)\in \partial \Omega
global  KM_V nod Ak xk nx ny N h1 h2
%% compute the Load Vector
gf = zeros((nx+1)*(ny+1),1);
S = h1*h2;
for k=1:N
        N1 = 1/4*(1-xk)*(1-xk)';
        N2 = 1/4*(1+xk)*(1-xk)';
        N3 = 1/4*(1+xk)*(1+xk)';
        N4 = 1/4*(1-xk)*(1+xk)';

        %%%the node index of this element 
        index_1 = nod(k,1);
        index_2 = nod(k,2);
        index_3 = nod(k,3);
        index_4 = nod(k,4);

        %%% approximation of u at the Guass point
        u_e = u(index_1)*N1+u(index_2)*N2+u(index_3)*N3+u(index_4)*N4;
        F_u = u_e.^3;

        Fe(1) = S/4*Ak*(F_u .*N1)*Ak';
        Fe(2) = S/4*Ak*(F_u .*N2)*Ak';
        Fe(3) = S/4*Ak*(F_u .*N3)*Ak';
        Fe(4) = S/4*Ak*(F_u .*N4)*Ak';
        for i=1:4
            L=nod(k,i);
            gf(L)=gf(L)+Fe(i);
        end
end
%% deal with the boundary condition
bi=[1:ny+1,1+(ny+1)*nx:(ny+1)*(nx+1),1+(ny+1)...
    :(ny+1):1+(ny+1)*(nx-1),2*(ny+1):(ny+1):(ny+1)*nx];
MM=2*(nx+ny);
for k=1:MM
    gf(bi(k))=0;
end

%% solve the result
psi=KM_V\gf;


            