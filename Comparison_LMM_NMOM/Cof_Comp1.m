clear all;
%%% Compare the efficiency of NMOM and LMM in computing ground state solutions of
%%% -\Delta u(x,y)  = g(x,y) |u(x,y)|^{p-1}u(x,y) , (x,y)\in \Omega = (-L,L)^2
%%% u(x,y) = 0, (x,y) \in \partial \Omega
cputime1 = cputime;
format long
%Args for model problem.
L =1 ;
global Kx Ky Nx Ny
a = -L; b = L;
c = -L; d = L;

%% discretization by sine pseudospectral method
Nx = 64; % mesh size
Ny = Nx;
global dx dy
dx = (b-a)/Nx; x = (a:dx:b)'; xj = x(2:Nx); %%only interior nodes
dy = (d-c)/Ny; y = (c:dy:d)'; yk = y(2:Ny); %%only interior nodes
[X,Y] = meshgrid(x,y); [Xj,Yk] = meshgrid(xj,yk);
global Tf
kx = pi/(b-a)*[1:Nx-1];
ky = pi/(d-c)*[1:Ny-1];
[Kx,Ky] = meshgrid(kx,ky); 
Tf=(Kx.^2+Ky.^2); 

%% initialization
U_0   = (1-Xj.^2).*(1-Yk.^2).*(2*(Xj-0.5).^2  + (Yk+0.5).^2);
U_0 = U_0/(inp(U_0,U_0)^(1/2));

%% fixed step-size
fix_step_temp = [1,0.1];

%% save the header
filename = ['eff_sensity_nmNMOM_nmLMM.xlsx'];
head = {'p', 'l', 'fix_step','iter_NMOM','CPU_NMOM', 'error_NMOM','energy_NMOM','iter_LMM', 'CPU_LMM', 'error_LMM','energy_LMM'}; 
xlswrite(filename,head);
D_index = 1;
epsilon = 1e-4;
    for p=1.5:0.5:3
        for l=p-1:0.5:p-0.5
            for j=1:length(fix_step_temp)
%           
                fix_step = fix_step_temp(j);  
                g_xy = (Xj.^2 +Yk.^2).^(l/2);
                
                tic
                [U_NMOM, energy_NMOM, iter_NMOM, error_temp_NMOM] = nm_NMOM(g_xy,U_0,fix_step,epsilon,p);
                cpu_time_NMOM = toc;
               
                tic
                [U_LMM, energy_LMM, iter_LMM, error_temp_LMM] = nm_LMM(g_xy,U_0,fix_step,epsilon,p);
                cpu_time_LMM = toc;
                D_index = D_index +1;
                
                xlswrite(filename,{p,l,fix_step,iter_NMOM,cpu_time_NMOM,error_temp_NMOM,energy_NMOM, iter_LMM,cpu_time_LMM, error_temp_LMM,energy_LMM},['A',num2str(D_index),':K',num2str(D_index)]);  

            end
        end
    end

%     

