function Plot_czx(u_1,energy_1)
%%% plot the profiles of the result
global nx ny x_left x_right y_down y_up h1 h2
U = zeros(nx+1,ny+1);
for i = 1:(nx+1)*(ny+1)
    col = mod(i,ny+1);
   if col ==0 
       col = ny+1;
   end
   row = (i-col)/(ny+1)+1;
   U(row,col)=u_1(i);
end

figure(1)
[X,Y]=meshgrid(x_left:h1:x_right,y_down:h2:y_up);
max_value = max(max(abs(U)));
[index_x,index_y] = find(U == max(max(U)));
max_x = X(index_x,index_y);
max_y = Y(index_x,index_y);

set(gcf,'Position',[0 150 300 280]);
set(gca,'Position',[.1,.1,.7,.75])
pcolor(X,Y,U);   shading interp;
title(['||u||_{\infty} = ',Get_deci(max_value,3) ,' at [',Get_deci(max_x,3) ,',',Get_deci(max_y,3),'], E(u) = ',Get_deci(energy_1,3)],'position',[0.1,1.2,0],'Fontsize',10,'FontWeight','Normal','FontName','Times New Roman');                                       % 标题需要更改
% my_title = get(gca,'title');
% set(my_title,'position',[0.1,1.2,0],'Fontsize',10,'FontWeight','Normal','FontName','Times New Roman')
set(gca,'FontSize',10,'FontWeight','Normal')
colormap('jet');                                     
hh = colorbar;
%POS = get(hh,'pos')
set(hh,'pos',[.85, .1 .05 .75]);