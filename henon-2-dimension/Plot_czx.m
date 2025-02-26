function  Plot_czx(u_1,energy_1)
global R Theta X Y
max_value = max(max(u_1));
[index_x,index_y] = find(u_1 == max_value);
if index_x == ones(size(index_x))
    max_x = 0;
    max_y = 0;
else
    max_x = R(index_x,index_y)*cos(Theta(index_x,index_y));
    max_y = R(index_x,index_y)*sin(Theta(index_x,index_y));
end
length_1 = length(R(1,:)); 
length_2 = length(R(:,1));
XX = zeros(length_2,length_1+1);
XX(:,1:length_1) = X;
XX(:,end) = X(:,1);
YY = zeros(length_2,length_1+1);
YY(:,1:length_1) = Y;
YY(:,end) = Y(:,1);
new_U = zeros(length_2,length_1+1);
new_U(:,1:length_1) = u_1;
new_U(:,end) = u_1(:,1);
figure(3);
set(gcf,'Position',[0 150 300 280]);
set(gca,'Position',[.1,.1,.7,.75])
pcolor(XX,YY,new_U);   shading interp;
title(['||u||_{\infty} = ',Get_deci(max_value,3) ,' at [',Get_deci(max_x,3) ,',',Get_deci(max_y,3),'], E(u) = ',Get_deci(energy_1,3)],'position',[0.1,1.2,0],'Fontsize',10,'FontWeight','Normal','FontName','Times New Roman');                                       % 标题需要更改
% my_title = get(gca,'title');
% set(my_title,'position',[0.1,1.2,0],'Fontsize',10,'FontWeight','Normal','FontName','Times New Roman')
set(gca,'FontSize',10,'FontWeight','Normal')
colormap('jet');                                     
hh = colorbar;
%POS = get(hh,'pos')
set(hh,'pos',[.85, .1 .05 .75]);
xlim([-1.1,1.1]); ylim([-1.1,1.1]);