function res = Sol_Poisson(f)
%%% solve the poisson equation -\Delta \psi(x,y) = f(x,y), (x,y)\in\Omega=\{(x,y): x^2+y^2<1\} 
%%% u(x,y) =0, (x,y) \in \partial \Omega
global M theta 
%% step 1: polar transform: x = rcos(\theta) + rsin (\theta) 
%%% transform the orignal equation to 
%%%-\psi_{rr} - 1/r \psi_r -1/r^2\psi_{\theta \theta} = f(r,\theta)

%% step 2 Fourier expansion in \theta: 
%%% \psi(r,\theta) = \sum_{m=0}^M(\psi_1m(r) cos(m\theta) + \psi_2m(r) sin(m\theta))
%%% f(r,\theta) = \sum_{m=0}^M(f_1m(r) cos(m\theta) + f_2m(r) sin(m\theta))
[f_1,f_2] = Comp_dfc(f); 

%% Step 3 transform the orignal equation to -\psi_1m''- 1/r \psi_1m' + (m/r)^2\psi_1m = f_1m, the same as \psi_2m
%%% Solve \psi_1m, \psi_2m 
global L_phi L_phi0 x_lgl w_lgl phi0 phi A B A0 T N
%%% Let v_1m(t) = \psi_1m( (t+1)/2 ), g_1m(t)= 1/4*(t+1)f((t+1)/2)
%%% X_N(m) = \{ v \in P_N, v(-1) = v(1) = 0\}
%%% basis function:  \phi_i = L_i - L_{i+2} £¨m\neq 0)
%%% basis function:  \phi_i = L_i - L_{i+1} £¨m =  0)
%%% weak formulation: ((t+1)v_1m', w') +m^2(1/(t+1) v_1m, w) = £¨I_Ng, w)

%%% m=0
g_10 = (x_lgl+1)/4.*f_1(:,1);  
g_20 = (x_lgl+1)/4.*f_2(:,1);
s_10=(T'*(g_10.*w_lgl)).*[[0:N-1]'+0.5;(N)/2];  %% g_10 = \sum_{j=0}^{N} s_10 (j) * L_j
s_20=(T'*(g_20.*w_lgl)).*[[0:N-1]'+0.5;(N)/2];  %% g_20 = \sum_{j=0}^{N} s_20 (j) * L_j¡£
F_1 =  L_phi0'*s_10; % (I_N g_10,\phi_i)
F_2 =  L_phi0'*s_20;   
x_1 = (A0)\F_1;
x_2 = (A0)\F_2 ;
psi_1(:,1) = phi0*x_1;
psi_2(:,1) = phi0*x_2;

%%% m \neq 0
for m=1:M
    g_1m = (x_lgl+1)/4.*f_1(:,m+1);  %% tr
    g_2m = (x_lgl+1)/4.*f_2(:,m+1);
    s_1m=(T'*(g_1m.*w_lgl)).*[[0:N-1]'+0.5;(N)/2];  %% g_1m = \sum_{j=0}^{N} s_1m (j) * L_j
    s_2m=(T'*(g_2m.*w_lgl)).*[[0:N-1]'+0.5;(N)/2];  %% g_2m = \sum_{j=0}^{N} s_2m (j) * L_j
    F_1 =  L_phi' *s_1m; % (I_N g_1m,\phi_i)
    F_2 = L_phi' * s_2m;  % (I_N g_2m,\phi_i)
    H = A+m^2*B;
    x_1 =H\F_1;
    x_2 = H\F_2;
    psi_1(:,m+1) = phi*x_1; 
    psi_2(:,m+1) = phi*x_2;
end
%% get the solution $\psi$
res = psi_1*(cos([0:M]'*(theta'))) +psi_2*(sin([0:M]'*(theta')));


