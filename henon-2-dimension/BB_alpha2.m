function [ alpha_new ] = BB_alpha2(u_0,u_1,Rgrad0,Rgrad1,k)
%%% compute the BB step-size
%%% if k is even  tau^1 = (s_k, s_k)/(s_k,v_k),
%%% else £¬tau^2 = (s_k,v_k)/(v_k,v_k);
s_k = u_1 -u_0;
v_k = Rgrad1 -Rgrad0;
if mod(k,2) ==0
    alpha_new = abs(inp( s_k,s_k)/inp( s_k,v_k));
else
    alpha_new = abs(inp( s_k,v_k)/inp( v_k,v_k));
end
