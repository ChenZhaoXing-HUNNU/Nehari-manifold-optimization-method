function [U_res, energy_res,iter,error] = LMM(g_xy,V_0,fix_step,epsilon,p)  
    %% compute the ground state solutions on Henon equation -\Delta u  = |x|^l |u|^{p-1}p by local minimax method
    global Tf dx dy
    %% fix the step-size
    step_k = fix_step;
%     rho = 0.85;
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
 
    %% compute the error
    L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
    error = max(max(abs(L_U_k)));
    
    %% enter the iterations
    k=1;
    
    while  error  > epsilon
        if k>100000
            break;
        end
        
        %% test the step-size 
        v = v_k - step_k*gradE_k;
        v_k =v/(inp(v,v)^(1/2)); %% update v_k
        up = inp(v_k,v_k);
        down = dx*dy*(sum(sum(g_xy.* abs(v_k.^(p+1) ))));
        t_k =  (up/down)^(1/(p-1));
        phi_k=t_k*v_k; %%compute the maximum point along v_k;
        energy_k = 0.5*inp(phi_k,phi_k) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_k.^(p+1) ))));
        varphi = idst2( dst2(g_xy.*abs(phi_k.^(p-1)).*phi_k)./(Tf)); %% 
        gradE_k  = phi_k - varphi;
        
        %% compute the error
        L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
        error = max(max(abs(L_U_k)));
        k=k+1;
    end
   U_res = phi_k;
   energy_res = energy_k;
   iter = k;
end