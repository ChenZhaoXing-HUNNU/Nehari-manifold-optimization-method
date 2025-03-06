function [U_res, energy_res,iter,error] = nm_NMOM(g_xy,U_0,fix_step,epsilon,p)  
    %%% compute the ground state solutions on Henon equation:
    %%% -\Delta u(x,y)  = g(x,y) |u(x,y)|^{p-1}u(x,y) , (x,y)\in \Omega = (-L,L)^2
    %%% u(x,y) = 0, (x,y) \in \partial \Omega
    %%% by Nehari Manifold Optimization Method
    %%% With nonmonotone step-size search rule
   global Tf dx dy
   %% initial  parameters for nonmonotone search
    sigma = 1e-3;
    beta = 0.25;
    rho = 0.85;
    
    %% initial guess
    up = inp(U_0,U_0);
    down = dx*dy*(sum(sum(g_xy.* abs(U_0.^(p+1) ))));
    phi_k= (up/down)^(1/(p-1))*U_0;
    energy_k = 0.5*inp(phi_k,phi_k) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_k.^(p+1) ))));

    %% compute the Riemannian gradient
    varphi = idst2( dst2(g_xy.*abs(phi_k.^(p-1)).*phi_k)./Tf);
    gradE_k  = phi_k - varphi; 
    gradG_k = 2*phi_k - (p+1)*varphi;
    R_gradE_k = gradE_k- (inp(gradE_k,gradG_k)/inp(gradG_k,gradG_k))*gradG_k;
    
     %% initial for nonmonoton rule
    Q_k = 1;
    C_k = energy_k;
    
    %% compute the error
    L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
    error = max(max(abs(L_U_k)));
    error_temp(1)  = error;
    k=1;
    tag_fail = 0;
    while  error  > epsilon
        if k>100000
            break;
        end
        step_k = fix_step;
        v = phi_k - step_k*R_gradE_k;
        up = inp(v,v);
        down = dx*dy*(sum(sum(g_xy.* abs(v.^(p+1) ))));
        phi_test= (up/down)^(1/(p-1))*v;
        energy_k_test = 0.5*inp(phi_test,phi_test) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_test.^(p+1) ))));
        while energy_k_test >1e-15+ C_k - sigma*step_k*inp(R_gradE_k,R_gradE_k)
            if step_k<1e-6
                tag_fail =1;
                break
            end
            step_k = step_k*beta;
            v = phi_k - step_k*R_gradE_k;
            up = inp(v,v);
            down = dx*dy*(sum(sum(g_xy.* abs(v.^(p+1) ))));
            phi_test= (up/down)^(1/(p-1))*v;
            energy_k_test = 0.5*inp(phi_test,phi_test) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_test.^(p+1) ))));
        end
        
        phi_k = phi_test;
        energy_k = energy_k_test;
       
        Q_k = rho*Q_k+1;
        C_k = C_k + (energy_k - C_k)/Q_k;
        
        %% compute the Riemannian gradient
        varphi = idst2( dst2(g_xy.*abs(phi_k.^(p-1)).*phi_k)./(Tf));
        gradE_k  = phi_k - varphi; 
        gradG_k = 2*phi_k - (p+1)*varphi;
        R_gradE_k = gradE_k- inp(gradE_k,gradG_k)/inp(gradG_k,gradG_k)*gradG_k;
        
        %% compute the error
        L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
        error = max(max(abs(L_U_k)));
        k=k+1;
    end
            %% return result
   U_res = phi_k;
   energy_res = energy_k;
   iter = k;
end
