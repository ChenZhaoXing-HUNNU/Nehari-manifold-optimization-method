function [K_BD,M_cBD] = Compu_K_M( )
%%% Compute the stiffness matrix and mass matrix 
%%% under linear finite element£»
%%% n: the numer of element
%%% h1£ºthe length of element
global n h
nod = zeros(n,2);
for j=1:n
        nod(j,:)=[j,j+1];  % node number
end
    
Ke =  1/h*[1,-1; -1,1]; %element stiffness matrix
Me = 1/6*h*[2 1;1 2]; %element mass matrix

ntriplets  = n*4; 
I_BD = zeros(ntriplets+2,1);
J_BD = zeros(ntriplets+2,1);
Mass_cBD = zeros(ntriplets+2,1); %% 
Stiff_BD = zeros(ntriplets+2,1); %% 

ntriplets = 0;
for k=1:n
      for i=1:2 
            for r=1:2
                ntriplets = ntriplets +1;
                I_BD(ntriplets) = nod(k,i);
                J_BD(ntriplets) = nod(k,r);
                if   I_BD(ntriplets) ==1 || I_BD(ntriplets) == n+1  
                        Mass_cBD(ntriplets) = 0;   
                        Stiff_BD(ntriplets) = 0;
               else
                    Mass_cBD(ntriplets) = Me(i,r);
                    if J_BD(ntriplets) ==1 || J_BD(ntriplets) == n+1
                        Stiff_BD(ntriplets) = 0;
                    else
                        Stiff_BD(ntriplets) = Ke(i,r);
                    end
                end
            end
       end
end

%% deal with the boundary condition£»
bi=[1,n+1];
for j =1:length(bi)
    ntriplets = ntriplets+1;
    I_BD(ntriplets) = bi(j);
    J_BD(ntriplets) = bi(j);
    Stiff_BD(ntriplets) = 1;
%     Mass_BD(ntriplets) = 1;
end
M_cBD = sparse (I_BD,J_BD,Mass_cBD,n+1,n+1); 
K_BD = sparse (I_BD,J_BD,Stiff_BD,n+1,n+1); 
end
