# NMOM  
This repository contains the Matlab codes used to reproduce the results presented in the paper 

Nehari manifold optimization and its application for finding unstable solutions of semilinear elliptic PDEs
 
Zhaoxing Chen, Wei Liu, Ziqing Xie, and Wenfan Yi 

[https://arxiv.org/abs/2404.09892](https://arxiv.org/abs/2404.09892)  

## Directory structure      
Each folder in this repository corresponds to a specific numerical test, and contains the scripts required to reproduce the results. Run the main matlab file in each folder, then the corresponding numerical results can be obtained directly.                   
- henon-1-dimension  
  + Contains codes for computing the ground state solution of  the Hénon equation in $\Omega = (-1,1) $ ,      
           
    $$
       \begin{cases}        
        u''(x) + |x|^{l} |u(x)|^{p-1}u(x) = 0,  &x \in \Omega,\\\                      
        u(x) = 0,   &x \in \partial \Omega.            
       \end{cases}                  
    $$
    
  + Files: 
    * Comp_henon_1d.m  --- the main matlab file   
    * Rie_grad.m --- compute the Riemannian gradient of the energy functional        
    * Retraction.m --- retraction mapping      
    * BB_alpha2.m --- compute the BB step-size  
    * Compu_KM --- get the stiffness matrix and mass matrix by FEM discretization        
 
- henon-2-dimension  
  + Contains codes for computing the ground state solution of the Hénon equation in $ \Omega =$ { $ (x,y) : x^2 + y^2 < 1 $ },                   
     
    $$
    \begin{cases}      
    \Delta u(x,y) + |x^2+y^2|^{l/2} |u(x,y)|^{p-1}u(x,y)  = 0,  &(x,y)\in \Omega, \\\                
     u(x,y) = 0,    &(x,y) \in \partial \Omega.   
    \end{cases}            
    $$
       
  + Files:
    * Comp_henon_1d.m  --- the main matlab file for compare nmNMOM and nmLMM
    * Comp_henon_2d.m  --- the main matlab file for compare NMOM and LMM
    * get_Amatrix.m --- get the required matrix by the spectral-Galerkin    
    * Rie_grad.m --- compute the Riemannian gradient of the energy functional        
    * Retraction.m --- retraction mapping    
    * LGL_pw.m --- Legendre-Gauss-Lobatto quadrature nodes and weights  
    * inp.m --- $H$-inner product $(\cdot, \cdot)_H $  
    * Sol_Poisson.m --- solve the poisson equation related to the Riemannian gradient  
    * Comp_dfc.m --- compute the Fourier coefficient      
    * Comp_dlc.m --- compute the coefficient under Legendre polynomial     
    * LegendreP.m ---compute the value n-Legendre in x    
    * Integ.m --- compute the integration in $ \Omega $      
    * BB_alpha2.m --- compute the BB step-size    
    * Plot_czx.m --- plot the profile of the solution    
    * Get_deci.m --- set the decimal place    

- NLSE-2-dimension
  + Contains the codes for computing the ground state solution of  the nonlinear Schrödinger equation in $\Omega = (-1,1)^2 $ ,      
    
    $$ 
    \begin{cases}      
     -\Delta u(x,y) + V(x,y)u(x,y)   = u^3(x,y) ,  &(x,y)\in \Omega, \\\              
     u(x,y)  = 0,   &(x,y) \in \partial \Omega.         
    \end{cases}
    $$
       
  + Files:        
    * Comp_NLSE_2d.m --- the main matlab file
    * get_KM.m --- get the stiffness matrix and mass matrix by FEM    
    * Rie_grad.m --- compute the Riemannian gradient of the energy functional      
    * Retraction.m --- retraction mapping      
    * inp.m --- $H$-inner product $(\cdot, \cdot)_H $  
    * Sol_Poisson.m --- solve the poisson equation related to the Riemannian gradient  
    * elestiff_V.m --- get the element mass matrix with variable coefficient $V(x,y)$      
    * gausspw.m --- get the Guass points and weights    
    * BB_alpha2.m --- compute the BB step-size          
    * Plot_czx.m --- plot the profile of the solution      
    * Get_deci.m --- set the decimal place  
    * V.m -- define the variable coefficient in NLSE     


- Comparison_LMM_NMOM    
  + Contains the codes for comparising the efficient of LMM and NMOM in computing the ground state solutions of Hénon equation in $\Omega = (-1,1)^2 $,            
  
    $$ 
    \begin{cases} \Delta u(x,y) + |x^2+y^2|^{l/2} |u(x,y)|^{p-1}u(x,y) = 0 ,   &(x,y)\in \Omega, \\\          
      u(x,y) = 0,   &(x,y) \in \partial \Omega.       
    \end{cases}               
    $$
    
  + Files:  
    * Cof_Comp2.m --- the main matlab file 
    * NMOM.m --- Nehari manifold optimization algorithm under fixed step-size and descent direction        
    * nm_NMOM.m --- Nehari manifold optimization algorithm under nonmonotone step-size and descent direction    
    * LMM.m --- local minimax under fixed step-size and descent direction      
    * nm_LMM.m --- local minimax algorithm under nonmonotone step-size and descent direction        
    * inp.m --- $H$-inner product $(\cdot, \cdot)_H $ 
    * dst2.m --- two-dimensional discrete sine transform                
    * idst2.m --- two-dimensional inverse discrete sine transform  

## Requirements
Please ensure that your version of MATLAB is R2016a or later before running the script
