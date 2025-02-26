function [phi,A,B,C,in_phi,T,L_phi] = get_Amatrix(m,N,t)
%%% phi(:,i) value of i-th basis function phi_i in Legendre-Guass-Lobatto 
%%% A(i,j) = \int_{-1}^1 phi_i' * phi_j dt
%%% B(i,j) = \int_{-1}^1 1/(t+1)* phi_i * phi_j dt
%%% C(i,j) = \int_{-1}^1 (t+1)* phi_i * phi_j dt
%%% in_phi(i,j) = \int_{-1}^1  phi_i * phi_j dt
%%%T(:,i) : value of i-lengendre polynomial ( L_i) in Legendre-Guass-Lobatto;
%%% L_phi(i,j)£º\int_{-1}^{1} L_i(t) * phi_j(t) dt, L_i : i-lengendre polynomial
for k=0:N
        T(:,k+1) = LegendreP(k,t); 
end
if m~=0
    for k =0:N-2
        phi(:,k+1) = LegendreP(k,t) - LegendreP(k+2,t);
    end
    L_phi = zeros(N+1,N-1);
        %% ¼ÆËã\int_I I_N g * \phi_k, k=1,2,\codts N-2
    for j =1:N-2
        L_phi(j,j) = 2/(2*(j-1)+1);
        L_phi(j+2,j) = -2/(2*(j-1)+5);
    end
    L_phi(N-1,N-1) = 2/(2*(N-2)+1);
    L_phi(N+1,N-1) = -2/N;
    A = zeros(N-1,N-1);
    B = zeros(N-1,N-1);
    C = zeros(N-1,N-1);
    in_phi = zeros(N-1,N-1);
    for i=0:N-2
        A(i+1,i+1) = 4*i+6;
        A(i+1,i+2) = 2*i+4;
        A(i+2,i+1) = 2*i+4;
        B(i+1,i+1) = 2*(2*i+3)/((i+1)*(i+2));
        B(i+1,i+2) = -2/(i+2);
        B(i+2,i+1) = -2/(i+2);
        C(i+1,i+4) = -2*(i+3)/((2*i+5)*(2*i+7));
        C(i+1,i+3) = -2/(2*i+5);
        C(i+1,i+2) = 2/((2*i+5)*(2*i+1))+2*(i+3)/((2*i+5)*(2*i+7));
        C(i+1,i+1) = 2/(2*i+1)+2/(2*i+5);
        C(i+3,i+1) = -2/(2*i+5);
        C(i+2,i+1) = 2/((2*i+5)*(2*i+1))+2*(i+3)/((2*i+5)*(2*i+7));
        C(i+4,i+1) =  -2*(i+3)/((2*i+5)*(2*i+7));
        in_phi(i+1,i+1) = 2/(2*i+1) +2/(2*i+5);
        in_phi(i+1,i+3) = -2/(2*i+5);
        if i-1 >0
            in_phi(i+1,i-1) = -2/(2*i+1);
        end
    end
    A(N:end,:) =[];
    A(:,N:end) = [];
    B(N:end,:) =[];
    B(:,N:end) = [];
    C(N:end,:) =[];
    C(:,N:end) = [];
    in_phi(:,N:end) = [];
    in_phi(N:end,:) = [];
elseif m==0
    for k =0:N-1
        phi(:,k+1) = LegendreP(k,t) - LegendreP(k+1,t);
    end
    L_phi = zeros(N+1,N);
        %% ¼ÆËã\int_I I_N g * \phi_k, k=1,2,\codts N-2
    for j =1:N-1
        L_phi(j,j) = 2/(2*(j-1)+1);
        L_phi(j+1,j) = -2/(2*(j-1)+3);
    end
    L_phi(N,N) = 2/(2*(N-1)+1);
    L_phi(N+1,N) = -2/N;
    A = zeros(N,N);
    C = zeros(N,N);
    in_phi = zeros(N,N);
    for i=0:N-1
    A(i+1,i+1) = 2*i+2;
    C(i+1,i+1) = 4*(i+1)/((2*i+1)*(2*i+3));
    C(i+1,i+2) = 4/((2*i+1)*(2*i+3)*(2*i+5));
    C(i+1,i+3) = -2*(i+2)/((2*i+3)*(2*i+5));
    C(i+2,i+1) = 4/((2*i+1)*(2*i+3)*(2*i+5));
    C(i+3,i+1) =-2*(i+2)/((2*i+3)*(2*i+5));
    in_phi(i+1,i+1) = 2/(2*i+1) +2/(2*i+3);
    in_phi(i+1,i+2) = -2/(2*i+3);
    if i>0
        in_phi(i+1,i) = -2/(2*i+1);
    end
    A(N+1:end,:) =[];
    A(:,N+1:end) = [];
    C(N+1:end,:) =[];
    C(:,N+1:end) = [];
    in_phi(N+1:end,:) =[];
    in_phi(:,N+1:end) = [];
    B = [];
    end
end