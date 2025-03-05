************************************************************************************************
This code compares the efficiency of NMOM and LMM in computing the ground state solution of  
        \Delta u(x,y) + |x^1+y^2|^{l/2} |u(x,y)|^{p-1}u(x,y) , (x,y)\in \Omega = (-L,L)^2
                                u(x,y) = 0,   (x,y) \in \partial \Omega
The main Matlab file is: Cof_Comp2.m 
Subroutines called are: 
                         NMOM.m----- Nehari manifold optimization algorithm under fixed step-size and descent direction
                         nm_NMOM.m --- Nehari manifold optimization algorithm under nonmonotone step-size and descent direction
                         LMM.m----- Local minimax under fixed step-size and descent direction
                         nm_LMM.m --- Local minimax algorithm under nonmonotone step-size and descent direction
                         inp.m --- Inner product ( , )_H
                         dst2.m --- Two-dimensional discrete sine Transform.
                         idst2.m --- Two-dimensional inverse discrete sine transform.