************************************************************************************************
This code computes the ground state solution of  the nonlinear Schr\"odinger equation : 
          - \Delta u(x,y) + V(x,y)u(x,y)   = u(x,y)^3 , (x,y)\in \Omega = (-1,1)^2
                                u(x,y) = 0,   (x,y) \in \partial \Omega

The main Matlab file is: Comp_NLSE_2d.m 
Subroutines called are: 
                        get_KM.m --- get the stiffness matrix and mass matrix by FEM
                        Rie_grad.m --- compute the Riemannian gradient of the energy functional
                        Retraction.m --- retraction mapping
                        inp.m --- the H-inner product.
                        Sol_Poisson.m --- solve the poisson equation -\Delta \psi(x,y) + V(x,y)\psi = u^3, (x,y)\in\Omega=\{(x,y): x^2+y^2<1\}
                                                 with Dirichlet boundary
                        elestiff_V.m --- get the element mass matrix with V(x,y)
                        gausspw.m --- get the Guass points and weights
                        BB_alpha2.m --- compute the BB step-size
                        Plot_czx.m --- plot the profile of the solution
                        Get_deci.m --- Set the decimal place
                        V.m -- define the variable coefficient in NLSE 

