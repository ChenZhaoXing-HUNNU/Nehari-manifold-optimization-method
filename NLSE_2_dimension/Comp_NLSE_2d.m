clear all
clc
%%%  compute the ground state solutions of  NLSE:
%-\Delta u(x,y)+(\omega*(x^2+y^2) +\gamma) u(x,y) = u^3(x,y) (x,y) \in\Omega = (-1,1)^2     
% (x,y) = 0, (x,y)\in \partial \Omega

%% discretization by FEM with bilinear rectangular element
global x_left x_right y_up y_down ny nx h1 h2 N
x_left = -1;
x_right = 1;
y_up = 1;
y_down = -1;
ny=100; %% the number of element at each row
nx=100; %% the number of element at each column
h1 = (x_right-x_left)/ny; % width of nodes
h2 = (y_up-y_down)/nx; % height of nodes
N=ny*nx;      % the number of  rectangle elements 

%%% Guass point and weights
global xk Ak
[xk,Ak]=gausspw(3);

%%% the node coordinates
x = zeros((ny+1)*(nx+1),1);
y = zeros((ny+1)*(nx+1),1);
for i=1:(nx+1)*(ny+1)
   col = mod(i,nx+1);  
   if col ==0 
       col = nx+1;
   end
   row = (i-col)/(nx+1)+1;
   x(i)  = x_left + (col-1)*h1 ;%% row
   y(i)  = y_down + (row-1)*h2 ;%% column
end

%%% mapping from element to node
global nod
for j=1:nx
    for i=(j-1)*ny+1:j*ny
        nod(i,:)=[i+j-1,i+j,i+j+ny+1,i+j+ny]; % the nodes
    end
end

%% parameters for NLSE
global omega lambda V_xy
omega = 45;
lambda = 4 ;
V_xy =  V(x,y);


%% stiffness matrix and mass matrix by FEM
global KM_V K 
[KM_V,K] = get_KM(); 
%%% KM_V: the sum of stiffness matrix and the mass matrix which contains V 
%%%  M: the stiffness matrix without V


%% initial point 
u_0 = (x-x_left).*(x-x_right).* (y-y_up).*(y-y_down);
[u_0,energy_0] = Retraction(0,0,u_0); % retraction map u_0 to N and returns the energy
RgradE0 = Rie_grad(u_0); % \nabla_{N} E(u_0), \nabla G(u_0)
norm_RgradE0 = (inp(RgradE0, RgradE0))^(0.5); % || \nabla E(u_0)||_H
%% parameters for nonmonotone step-size search rule
Q_1 = 1;
alpha_max =10;
alpha_min=1;
C_1 = energy_0; 
beita = 0.25;
sigma = 0.001; 
varrho = 0.85;
t_0 = 1;


%% start iteration
tag = 0; % tag for the failure in finding step-size   
[u_k,energy_k] = Retraction(t_0, -RgradE0,u_0);
 while energy_k >  1e-15 + C_1 - sigma*t_0* inp( RgradE0,RgradE0) 
    t_0 = t_0*beita;
    [u_k,energy_k] = Retraction(t_0, -Rgrad0, u_0);
    if t_k<1e-4
        tag=1; 
        break;
    end
 end
 if tag ==1
    error('can not fund the step-size when k=1')
end

u_1 = u_k;
energy_1 = energy_k;
RgradE1= Rie_grad(u_1);
norm_RgradE1 = (inp(RgradE1,RgradE1))^(0.5);
Q_1 = varrho * Q_1 +1;
C_1 = ((Q_1-1)*C_1 + energy_1)/Q_1;

tag =0; %tag for the failure in finding step-size  
k=2;
while norm_RgradE1 >1e-6
    alpha_1= BB_alpha2(u_0, u_1, RgradE0,RgradE1,k); % BB step-size
   
    if alpha_1 <alpha_min
        alpha_1=alpha_min;
    elseif alpha_1>alpha_max
        alpha_1=alpha_max;
    end
    
    [u_k,energy_k] = Retraction(alpha_1, -RgradE1,u_1);
     t_k = alpha_1;
     while energy_k >  1e-15 + C_1 - sigma*t_k* inp( RgradE1,RgradE1) 
         if t_k<1e-5
            tag=1;
            break;
         else
            t_k = t_k*beita;
            [u_k,energy_k] = Retraction(t_k, -RgradE1, u_1);
         end 
     end
    
    %%% update information
    u_0= u_1;
    energy_0 = energy_1;
    RgradE0 = RgradE1; 
    norm_RgradE0 = norm_RgradE1; 
    
    u_1 = u_k;
    energy_1 = energy_k;
    RgradE1= Rie_grad(u_1);
    norm_RgradE1 = (inp(RgradE1,RgradE1))^(0.5);
    
    Q_1 = varrho * Q_1 +1;
    C_1 = ((Q_1-1)*C_1 + energy_1)/Q_1;
    k=k+1;
    
    if (k>500) || (tag==1)
        break;
    end
end

u_res =u_1;
energy_res = energy_1;

Plot_czx(u_1,energy_1)
