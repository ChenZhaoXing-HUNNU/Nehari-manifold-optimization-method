function inp_uv = inp(u,v)
%%%compute the inner product \ing_{Omega} \nabla u \nabla v dxdy 
%%% = \int_0^{2pi} \int_{0}^1 r u_r*v_r dr d \theta +\int_0^1 \int_{0}^{2pi} 1/r u_{\theta}v_{\theta} drd \theta
global   A  B  A0
global M 

%% compute the Fourier cofficient of u,v
[u_1m,u_2m] = Comp_dfc(u); 
[v_1m,v_2m] = Comp_dfc(v); 

%% compute the coefficient of u_1m expanded on X_u = {\phi_0,\phi_1,\cdots\}
[u_X_10,u_X_1] = Comp_dlc(u_1m);
[u_X_20,u_X_2] = Comp_dlc(u_2m);
[v_X_10,v_X_1] = Comp_dlc(v_1m);
[v_X_20,v_X_2] = Comp_dlc(v_2m);
inp_uv = 0;
for m=0:M
    if m==0
        inp_uv = inp_uv+2*pi*(u_X_10'*A0*v_X_10+u_X_20'*A0*v_X_20);
    else    
        inp_uv = inp_uv+pi*((u_X_1(:,m)'*A*v_X_1(:,m) +u_X_2(:,m)'*A*v_X_2(:,m))+m^2*(u_X_1(:,m)'*B*v_X_1(:,m)+u_X_2(:,m)'*B*v_X_2(:,m)));
    end
end
end


