function [ALL_left_BD,K_BD] = get_KM()
%%% compute the stiffness matrix and mass matrix
%%% ALL_left_BD: the sum of stiffness matrix and the mass matrix with V. 
%%% K_BD : the stiffness matrix
global  h1 h2 nod nx ny

Ke = (1/6)*((h2/h1)*[2 -2 -1 1;-2 2 1 -1;-1 1 2 -2;1 -1 -2 2]+(h1/h2)*[2 1 -1 -2;1 2 -2 -1;-1 -2 2 1;-2 -1 1 2]);
ntriplets  = ny*nx*16; % number of elements
I_BD = zeros(ntriplets+2*(nx+ny),1);
J_BD = zeros(ntriplets+2*(nx+ny),1);
K = zeros(ntriplets+2*(nx+ny),1); % Mass matrix without V(x,y)
ALL_BD = zeros(ntriplets+2*(nx+ny),1); % the sum of the 'Mass matrix with V(x,y)' and the 'stiffness matrix'
ntriplets = 0;
for k=1:nx*ny
      Me_V = elestiff_V(k);
      for i=1:4 
            for r=1:4 
                ntriplets = ntriplets +1;
                I_BD(ntriplets) = nod(k,i);
                J_BD(ntriplets) = nod(k,r);
                if mod(I_BD(ntriplets), nx+1)==0 || mod(I_BD(ntriplets)-1, nx+1)==0  || I_BD(ntriplets)<nx+2 || I_BD(ntriplets)>(nx+1)*ny 
                    K(ntriplets) = 0; % deal with the Dicrichlet boundary£»
                    ALL_BD(ntriplets) = 0; % deal with the Dicrichlet boundary£»
                else

                    if mod(J_BD(ntriplets), nx+1) ==0 || mod(J_BD(ntriplets)-1, nx+1) ==0 || J_BD(ntriplets)<nx+2 || J_BD(ntriplets)>(nx+1)*ny
                        ALL_BD(ntriplets) = 0; % deal with the Dicrichlet boundary£»
                        K(ntriplets) = 0; 
                    else 
                        ALL_BD(ntriplets) = Ke(i,r)+ Me_V(i,r);
                        K(ntriplets) = Ke(i,r);
                    end
                end
            end
      end
end

%%%  deal with the Dicrichlet boundary£»

%% deal with the boundary condition
bi=[1:nx+1,1+(nx+1)*ny:(nx+1)*(ny+1),1+(nx+1)...
    :(nx+1):1+(nx+1)*(ny-1),2*(nx+1):(nx+1):(nx+1)*ny];
for j =1:length(bi)
    ntriplets = ntriplets+1;
    I_BD(ntriplets) = bi(j);
    J_BD(ntriplets) = bi(j);
    ALL_BD(ntriplets) = 1;
    K(ntriplets) = 1;
end
K_BD = sparse (I_BD,J_BD,K,(nx+1)*(ny+1),(nx+1)*(ny+1)); 
ALL_left_BD = sparse (I_BD,J_BD,ALL_BD,(nx+1)*(ny+1),(nx+1)*(ny+1)); 