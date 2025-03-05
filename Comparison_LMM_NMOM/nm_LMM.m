function [U_res, energy_res,iter,error] = nm_LMM(g_xy,V_0,fix_step,epsilon,p)  
    %%% compute the ground state solutions on Henon equation:
    %%% -\Delta u(x,y)  = g(x,y) |u(x,y)|^{p-1}u(x,y) , (x,y)\in \Omega = (-L,L)^2
    %%% u(x,y) = 0, (x,y) \in \partial \Omega
    %%% by local minimax method
    %%% With nonmonotone step-size search rule
    global Tf dx dy
    %%  initial  parameters for nonmonotone search
    rho = 0.85;
    sigma = 1e-3;
    beta = 0.25;
    
   %% initial set U_0 \in S
    v_k = V_0/(inp(V_0,V_0)^(1/2));
    up = inp(v_k,v_k);
    down = dx*dy*(sum(sum(g_xy.* abs(v_k.^(p+1) ))));
    t_k =  (up/down)^(1/(p-1));
    phi_k  = t_k*v_k; %% compute the maximum point 
    energy_k = 0.5*inp(phi_k,phi_k) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_k.^(p+1) ))));
   
    %% compute the H-gradient at phi_k
    varphi = idst2( dst2(g_xy.*abs(phi_k.^(p-1)).*phi_k)./(Tf)); %% 
    gradE_k  = phi_k - varphi;
 
    %% initial for nonmonoton rule
    Q_k = 1;
    C_k = energy_k;
    
    %% compute the error
    L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
    error = max(max(abs(L_U_k)));
    error_temp(1)  = error;
   
    %% enter the iterations
    k=1;
    tag_fail = 0;
  
    while  error  > epsilon
        if k>100000
            break;
        end
        step_k = fix_step;
        v = v_k - step_k*gradE_k;
        v_test =v/(inp(v,v)^(1/2)); %% update v_k
        
        up = inp(v_test,v_test);
        down = dx*dy*(sum(sum(g_xy.* abs(v_test.^(p+1) ))));
        t_test =  (up/down)^(1/(p-1));
        phi_test=t_test*v_test; %%compute the maximum point along v_k;
        energy_k_test = 0.5*inp(phi_test,phi_test) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_test.^(p+1) ))));
        while energy_k_test > 1e-15 + C_k -sigma*step_k*t_k*inp(gradE_k,gradE_k)
             if step_k <1e-6
                 tag_fail =1;
                 break
             end
             step_k = step_k*beta;
             v = v_k - step_k*gradE_k;
             v_test =v/(inp(v,v)^(1/2)); %% update v_k
             up = inp(v_test,v_test);
             down = dx*dy*(sum(sum(g_xy.* abs(v_test.^(p+1) ))));
             t_test = (up/down)^(1/(p-1));
             phi_test= t_test*v_test; %%compute the maximum point along v_k;
             energy_k_test = 0.5*inp(phi_test,phi_test) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_test.^(p+1) ))));
        end
        if tag_fail ==1
            fprintf('LMM can not find the step-size\n');
            break
        end
        v_k = v_test;
        phi_k = phi_test;
        t_k = t_test;
        energy_k = energy_k_test;
        
        varphi = idst2( dst2(g_xy.*abs(phi_k.^(p-1)).*phi_k)./(Tf)); %% 
        gradE_k  = phi_k - varphi;
        
        Q_k = rho*Q_k+1;
        C_k = C_k + (energy_k - C_k)/Q_k;
        %% compute the error
        L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
        error = max(max(abs(L_U_k)));
        k=k+1;

    end
   U_res = phi_k;
   energy_res = energy_k;
    iter = k;
end