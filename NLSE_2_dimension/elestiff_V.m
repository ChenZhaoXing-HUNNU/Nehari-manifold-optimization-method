function [Me_V ]= elestiff_V( k )
%%% compute the element mass matrix with V(x,y)
%%% use Guass integration \int_{e_i) V*N*N'dx
global h1 h2 xk Ak ny x_left y_down
S = h1*h2;   
ele_id = k;
col = mod(ele_id,ny); %%判断单元处于哪一列
    if col==0
        col = ny;
    end
row = (ele_id-col)/ny+1; 
cx = col-1;                    
dx = row-1; 

%%% compute the value of V at Guass point。
x_c = x_left + cx * h1 +0.5*h1; 
y_c = y_down + dx * h2 +0.5*h2;
x = (h1/2*xk+x_c)*ones(1,length(xk)); 
y = ones(length(xk),1)*(h2/2*xk+y_c)';
V_ek = V(x,y);

%%% the value of basis function at Guass point。
N = zeros(length(xk),length(xk),4);
N(:,:,1) = 1/4*(1-xk)*(1-xk)';
N(:,:,2) = 1/4*(1+xk)*(1-xk)';
N(:,:,3) = 1/4*(1+xk)*(1+xk)';
N(:,:,4) = 1/4*(1-xk)*(1+xk)';


%%
Me_V = zeros(4,4);
for ii=1:4
    for jj = 1:4
        Me_V(ii,jj) = Ak*(N(:,:,ii).*N(:,:,jj).*V_ek)*Ak';
    end
end
Me_V = S/4*Me_V;


