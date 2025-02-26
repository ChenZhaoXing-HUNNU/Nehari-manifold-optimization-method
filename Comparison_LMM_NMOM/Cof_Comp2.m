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
fix_step_temp = [1, 0.1, 0.01, 0.001];
D_index = 1;
epsilon = 1e-4;
    for p=1.5
        for l=0.5
            for j=1:length(fix_step_temp)
%              for j=1:4
                fix_step = fix_step_temp(j);  
                g_xy = (Xj.^2 +Yk.^2).^(l/2);
                tic
%                 [U_NMOM, energy_NMOM, iter_NMOM, error_temp_NMOM] = NMOM(g_xy,U_0,fix_step,epsilon,p);
%                 cpu_time_NMOM = toc;
%                 
%                   tic
                [U_NMOM, energy_NMOM, iter_NMOM, error_temp_NMOM] = nm_NMOM(g_xy,U_0,fix_step,epsilon,p);
                cpu_time_NMOM = toc;
                
%                  xlswrite(filename,{p,l, iter_NMOM,cpu_time_NMOM,error_temp_NMOM(end),energy_NMOM, iter_LMM,cpu_time_LMM, error_temp_LMM(end),energy_LMM},['A',num2str(D_index),':J',num2str(D_index)]);  

            % %     save(['eff_comp//test//m-',num2str(m),'//c-1//alpha_0_01//error_m-',num2str(m),'_RNAG_N_',num2str(beta),'.mat'],'error_temp_nag')
            %     D_index = D_index+1;
                        % %     
%                 tic
%                 [U_LMM, energy_LMM, iter_LMM, error_temp_LMM] = LMM(g_xy,U_0,fix_step,epsilon,p);
%                 cpu_time_LMM = toc;
                
%                  tic
                [U_LMM, energy_LMM, iter_LMM, error_temp_LMM] = nm_LMM(g_xy,U_0,fix_step,epsilon,p);
                cpu_time_LMM = toc;
                D_index = D_index +1;
%                 xlswrite(filename,{p,l,fix_step,iter_NMOM,cpu_time_NMOM,error_temp_NMOM,energy_NMOM, iter_LMM,cpu_time_LMM, error_temp_LMM,energy_LMM},['A',num2str(D_index),':K',num2str(D_index)]);  

        %     save(['eff_comp//test//m-',num2str(m),'//c-1//alpha_0_01//error_m-',num2str(m),'_MRNAG_N_',num2str(beta),'.mat'],'error_temp_nag_nonm')
        %     D_index = D_index+1;
        %     xlswrite(filename,{'RAG_nonm_N',beta, fix_step,iter_nag_nonm, cpu_time_nag_nonm,energy_nag_nonm, error_temp_nag_nonm(end)},['A',num2str(D_index),':G',num2str(D_index)]);  
                %     
            end
        end
    end

%     

