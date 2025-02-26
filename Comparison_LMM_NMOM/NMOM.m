function [U_res, energy_res,iter,error] = NMOM(g_xy,U_0,fix_step,epsilon,p)  
    %% compute the ground state solutions on Henon equation -\Delta u  = |x|^l |u|^{p-1}p
    global Tf dx dy
   %% fix the step-size
    step_k = fix_step;

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
     
    %% compute the error
    L_U_k = idst2(Tf.*dst2(phi_k))  - g_xy.*abs(phi_k.^(p-1)).*phi_k;
    error = max(max(abs(L_U_k)));
    k=1;

    while  error  > epsilon
        if k>100000
            break;
        end
        v = phi_k - step_k*R_gradE_k;
        up = inp(v,v);
        down = dx*dy*(sum(sum(g_xy.* abs(v.^(p+1) ))));
        phi_k= (up/down)^(1/(p-1))*v;
        energy_k= 0.5*inp(phi_k,phi_k) - 1/(p+1)*dx*dy*(sum(sum(g_xy.* abs(phi_k.^(p+1) ))));
       
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
            
   U_res = phi_k;
   energy_res = energy_k;
   iter = k;
end