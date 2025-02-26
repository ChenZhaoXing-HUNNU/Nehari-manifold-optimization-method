function [bp,wf]=gausspw(n) 
%A standard  code to get the gauss points and the weights
    u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
    [vc,bp]=eig(diag(u,-1)+diag(u,1));
    [bp,k]=sort(diag(bp));% the integration points
    wf=2*vc(1,k).^2;%The weight