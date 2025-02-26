function integ_f = integ(f)
%%% compute \int_{\Omega} f dxdy
global w_lgl R
tran_f = f.*R; % dxdy = r drd\theta -- R Jacobi matrix
f_theta = 2*pi/length(f(1,:))*(sum(tran_f,2));  % \theta -coordinate
integ_f = 0.5*w_lgl'*f_theta; % r-coordinate , r = (t+1)/2