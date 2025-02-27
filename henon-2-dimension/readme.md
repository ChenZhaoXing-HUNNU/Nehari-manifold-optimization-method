************************************************************************************************
This code computes the ground state solution of  the H\'enon equation : 
          \Delta u(x,y) + |x^2+y^2|^{l/2} |u(x,y)|^{p-1}u(x,y) = 0, (x,y)\in \Omega = \{(x,y):x^2+y^2<1\}   
                                u(x,y) = 0,   (x,y) \in \partial  

The main Matlab file is: Comp_henon_2d.m   
Subroutines called are:   
                        get_Amatrix.m --- get the required matrix by the spectral-Galerkin  
                        Rie_grad.m --- compute the Riemannian gradient of the energy functional  
                        Retraction.m --- retraction mapping
                        LGL_pw.m --- Legendre-Gauss-Lobatto quadrature nodes and weights
                        inp.m --- the H-inner product.
                        Sol_Poisson.m --- solve the poisson equation -\Delta \psi(x,y) = f(x,y), (x,y)\in\Omega=\{(x,y): x^2+y^2<1\}
                                                 with Dirichlet boundary
                        Comp_dfc.m --- compute the Fourier coefficient
                        Comp_dlc.m --- compute the coefficient under Legendre polynomial 
                        LegendreP.m ---compute the value n-Legendre in x
                        Integ.m --- compute the integration  
                        BB_alpha2.m --- compute the BB step-size
                        Plot_czx.m --- plot the profile of the solution
                        Get_deci.m --- Set the decimal place

