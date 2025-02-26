function [X_0, X] = Comp_dlc(u)
%%% compute the coefficient of u = \{u_m\}_{m=1}^M under the basis \{\phi\}
%%% when m=0 £¬ basis function \phi_j = L_j - L_{j+1}, j=0,1,\cdots,N-1. 
%%% when m\neq 0£¬basis function \phi_j = L_j -L_{j+2}, j=0,1,\dots, N-2
global M w_lgl in_phi0  in_phi L_phi L_phi0 T N

for m=0:M
    if m==0
        s_m=(T'*(u(:,m+1).*w_lgl)).*[[0:N-1]'+0.5;(N)/2];  %% compute the coefficient of u_m under Legendre basis \{L_j\}
        now_f_1 =  L_phi0'*s_m;  % compute (I_Nu_m, \phi_k)
        X_0 = (in_phi0)\now_f_1; 
    else
        s_m=(T'*(u(:,m+1).*w_lgl)).*[[0:N-1]'+0.5;(N)/2];
        now_f_1 =  L_phi'*s_m; % compute (I_Nu_m, \phi_k)
        X(:,m) = in_phi \now_f_1;
    end
end