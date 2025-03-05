clear all;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
format long;
%%% Use NMOM to compute the ground state solution of H\'enon equation: 
%%% \Delta u(x,y) + sqrt{x^2+y^2}^{ l } |u(x,y)|^{p-1}u(x,y)=0, (x,y) \in \Omega = {(x,y):x^2+y^2<1}
%%% u(x,y) = 0; (x,y) \in \partial \Omega
%% Discretization of SPECTRAL-GALERKIN METHODS
global M N
M =64 ; % u(r,\theta) = sum_{m=0}^{N}(u^{1m}(r)cos(m\theta)+u^{2m}(r)sim(m\theta));
N = 64 ; % the max degree of polynomials for u^{1m}(r), u^{2m}(r), basis function phi_i
global x_lgl w_lgl
[x_lgl,w_lgl] = LGL_pw(N); % N-Legendre-Guass-Lobatto points and weights
global r theta
r = (x_lgl+1)/2; 
theta = (pi*[0:2*M-1]/M)';  % polar transformation 
global Theta R X Y
[Theta,R] = meshgrid(theta,r);
X = R.*cos(Theta);
Y = R.*sin(Theta);
global phi in_phi A  B C T L_phi 
%%% phi(:,i) value of i-th basis function phi_i in Legendre-Guass-Lobatto 
%%% A(i,j) = \int_{-1}^1 phi_i' * phi_j dt
%%% B(i,j) = \int_{-1}^1 1/(t+1)* phi_i * phi_j dt
%%% C(i,j) = \int_{-1}^1 (t+1)* phi_i * phi_j dt
%%%T(:,i) : value of i-lengendre polynomial ( L_i) in Legendre-Guass-Lobatto;
%%% L_phi(i,j)£º\int_{-1}^{1} L_i(t) * phi_j(t) dt, L_i : i-lengendre polynomial
[phi,A,B,C,in_phi,T,L_phi] = get_Amatrix(1,N,x_lgl) ; % phi_i =L_i - L_{i+2} 
global phi0 in_phi0 A0  B0 C0 T_0 L_phi0
[phi0,A0,B0,C0,in_phi0,T_0, L_phi0] = get_Amatrix(0,N,x_lgl); % phi_i = L_i - L_i{i+1}

%% parameters for H\'enon equation
global p l g_l
p=2 ;
l= 2;
g_l =  R.^l;

%% initial point in N for NMOM
u_0= -(R.^2-1).*exp(2*R.*cos(Theta)+R.*sin(Theta)); 
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

u_result = u_1;
energy_result = energy_1;
%% plot the profile of the result u
Plot_czx(u_result,energy_result)

