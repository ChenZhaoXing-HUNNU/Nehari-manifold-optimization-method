************************************************************************************************
This code computes the ground state solution of  the H\'enon equation : 
         u''(x) + |x|^{l} |u(x)|^{p-1}u(x) = 0, (x)\in \Omega = (-1,1)
                                u(x) = 0,   (x) \in \partial \Omega
************************************************************************************************
The main Matlab file is: Comp_henon_1d.m 
Subroutines called are: 
                        Rie_grad.m --- compute the Riemannian gradient of the energy functional
                        Retraction.m --- retraction mapping
                        BB_alpha2.m --- compute the BB step-size
                        Compu_KM --- get the stiffness matrix and mass matrix by FEM discretization
                     

