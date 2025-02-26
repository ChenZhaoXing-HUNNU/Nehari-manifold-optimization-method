function [ alpha_new ] = BB_alpha2( u_0,laplace_u_0,u_1,laplace_u_1,gradE0,laplace_gradE0,gradE1,laplace_gradE1,k)
%%% compute the BB step-size
%%% if k is even  tau^1 = (s_k, s_k)/(s_k,v_k),
%%% else £¬tau^2 = (s_k,v_k)/(v_k,v_k);
s_k = u_0 - u_1;
laplace_sk  = laplace_u_0 - laplace_u_1;
v_k = gradE0 - gradE1;
laplace_vk = laplace_gradE0 - laplace_gradE1;
if mod(k,2)==0
    alpha_new =  abs((laplace_sk'*s_k )/(laplace_sk'*v_k));  
else
    alpha_new =  abs((laplace_sk'*v_k )/(laplace_vk'*v_k));  
end
