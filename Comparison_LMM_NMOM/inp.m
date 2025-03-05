function inp_f_g = inp(f,g) 
%%% compute the inner prodcut (f,g)_H.
    global Tf dx dy
        inp_f_g = dx*dy*sum(sum(idst2(Tf.*dst2(f)).*g));

end