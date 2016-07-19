1. Test files: 
      a. test_lorenz.m
      b. test_rossler.m
      c. test_hyperchaos.m

2. Algorithm: am_solver.m  
     Alternating minimization to solve the following joint sparsity + sparsity problem
             (C,E) =  min mu/2||phiX*C + E - Xdot + b||_F^2 +  sum_j||E^j||_2
                    subject to C is sparse
                b  =  b + phiX*C + E - Xdot 


3. Supplementary Material:
      a. corrupted_data.m:   add bandwith corruption to clean data X
      b. dictionary3.m:      construct dictionary matrix built from the 3D data
      c. dictionary4.m:      construct dictionary matrix built from the 4D data
      c. time_derivative.m:  1st/2nd order approximation of time derivative 

%====================================================================================
Authors: Giang Tran and Rachel Ward
Institution: The University of Texas at Austin
Version: 1.0, Year: 2016
More information can be found at: 
    G. Tran and R. Ward, "Exact Recovery of Chaotic Sysmtems from Highly Corrupted Data", 
                         https://arxiv.org/abs/1607.01067 
