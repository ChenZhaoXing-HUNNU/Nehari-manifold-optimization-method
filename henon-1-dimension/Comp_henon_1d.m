clear all;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
format long;
%%% Use NMOM to compute the ground state solution H\'enon equation:
%%% u'' + |x|^{ l }|u|^{ p-1 } u =0, x \in (-L, L)
%%%  u(-L） = u(L) =0

%% The domain $\Omega = (-L, L)$
global x_left x_right
L = 1; 
x_left = -L;
x_right = L;

%% Discretization of FEM
global n h x
n = 1000;    % the number of elements
h = (x_right-x_left)/n;   % the length of element
x = [x_left:h:x_right]'; % node coordinates 

%% paramters in H\'enon equation 
global p l g_xy
p = 2;
l = 4;
g_xy = abs(x).^l;
 
%% Prepare the stiffness matrix and the mass matrix after discretization
%%% M_cBD：mass matrix，
%%% K_BD ：Stiffness matrix
global K_BD M_cBD
[K_BD,M_cBD] = Compu_K_M();

%% initial point u_0 and \nabla_{N} E(u_0) 
u_0 = (x-x_left).*((x-x_right).^2);
laplace_u0 = 6*x-2;
[u_0,laplace_u0,energy_0] = Retraction(0, 0,u_0, laplace_u0, 0); 
[RgradE0, laplace_RgradE0]= Rie_grad(u_0, laplace_u0); % \nabla_{N} E(u_0)
norm_RgradE0 = sqrt(-h*laplace_RgradE0'*RgradE0); %\| \nabla_{N} E(u_0)\|_H

%% start use NMOM minimizing $E$ on N.
%%% initialize the parameters in NMOM
Q_1 = 1;
C_1 = energy_0; 
beita = 0.25;
sigma = 0.001; 
varrho = 0.85;
epsilon = 1e-6; 

%%
%%% 
t_0 =1; %initial step-size
[u_k,laplace_uk,energy_k] = Retraction(t_0,RgradE0,u_0,laplace_u0,laplace_RgradE0);
tag = 0; % tag for the failure in finding step-size   
while energy_k >  1e-6 + C_1 - sigma*t_0* (-h*laplace_RgradE0'*RgradE0)
     t_0 = t_0 *beita;
     [u_k,laplace_uk,energy_k] = Retraction(t_0,RgradE0,u_0,laplace_u0,laplace_RgradE0);
     if t_0 <1e-4
         tag =1;
         break;
     end
end
if tag ==1
    error('can not fund the step-size when k=1')
end

u_1 = u_k;
laplace_u1 = laplace_uk;
energy_1 = energy_k;
[RgradE1, laplace_RgradE1]= Rie_grad(u_1, laplace_u1);
norm_RgradE1 = sqrt(-h*laplace_RgradE1'*RgradE1);
Q_1 = varrho * Q_1 +1;
C_1 = ((Q_1-1)*C_1 + energy_1)/Q_1;

tag =0; %tag for the failure in finding step-size
k =2;
while (norm_RgradE1 >epsilon)
%     alpha_1 = BB_alpha( h1,h2,t_0,gradE0,gradE1,u_0,g_xy,V_xy,K,p)   
    alpha_1 = BB_alpha2( u_0,laplace_u0,u_1,laplace_u1,RgradE0,laplace_RgradE0,RgradE1,laplace_RgradE1,k);
    if alpha_1 < 1  
        alpha_1 = 1;
    end
    if alpha_1 >10
            alpha_1 = 10;

    end
 
    [u_k,laplace_uk,energy_k] = Retraction(alpha_1,RgradE1,u_1,laplace_u1,laplace_RgradE1);%%由拉回算子计算得到第一步。
    t_k = alpha_1;
    
    while energy_k >  1e-15 + C_1 - sigma*t_k* (-h*laplace_RgradE1'*RgradE1) 
        if t_k <1e-5
            tag=1;            
            break
         else
             t_k = t_k *beita;
             [u_k,laplace_uk,energy_k] = Retraction(t_k,RgradE1,u_1,laplace_u1,laplace_RgradE1);
         end
    end
    %%
    if tag==1 || k>5000
        break;
    end
    k = k+1;
    
    %%% update 
    
    u_0 = u_1; 
    laplace_u0 = laplace_u1;
    energy_0 = energy_1;
    RgradE0 = RgradE1; 
    laplace_RgradE0 = laplace_RgradE1;
    
    u_1 = u_k;
    laplace_u1 = laplace_uk;
    energy_1 = energy_k;
    [RgradE1, laplace_RgradE1]= Rie_grad(u_1, laplace_u1);
    norm_RgradE1 = sqrt(-h*laplace_RgradE1'*RgradE1);
   
    Q_1 = varrho * Q_1 +1;
    C_1 = ((Q_1-1)*C_1 + energy_1)/Q_1;
end

%% obtaine the final ground state solution
result_u = u_1;
laplace_result_u = laplace_u1;

