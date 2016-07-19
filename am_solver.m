%==========================================================================
% am_solver: alternating minimization to solve the following joint sparsity + sparsity problem
%             (C,E) =  min mu/2||phiX*C + E - Xdot + b||_F^2 +  sum_j||E^j||_2
%                    subject to C is sparse
%                b  =  b + phiX*C + E - Xdot 
%
% Input: 
%       Xdot(m x n): time derivative of data, kth row is the measurement value at time k*dt 
%       phiX(m x r): dictionary matrix built from the data 
%                         where n = dimension of the ODE system
%                               r = # basis functions
%                               m = # measurements
%       opts.      : An optional struct of options
%            tol   : tolerance for relative change of the outliers, default = 5e-3
%            maxit : max number of iterations, default = 100
%            coef_thres: thresholding parameter for the coefficient (coeffs - hard thresholding)
%            mu    : weight of the term ||phiX*C + E - Xdot + b||_F^2
% Output: 
%       coeff      : Recovered coefficients C
%       outlier    : Values of outliers E
%       opts.      : A complete struct of options, containing all the values 
%                        that were used by the solver.                        
%       outs.      : A struct with convergence information
%            iter  : #iterations
%            relerror_coeff: norm2 of two consecutive coeff's iterates each iteration
%            relerror_outlier: norm2 of two consecutive outlier's iterates each iteration
%            len_outlier: #nonzero entries in one column of outlier E each iteration

%
% Note: #rows(Xdot) = #rows(phiX), but #rows(Xdot) < #rows(X)
%
% Authors: Giang Tran and Rachel Ward
% Institution: The University of Texas at Austin
% Version: 1.0, Year: 2016
% More information can be found at: 
%     G. Tran and R. Ward, "Exact Recovery of Chaotic Sysmtems from Highly
%     Corrupted Data", https://arxiv.org/abs/1607.01067 
%==========================================================================

function [coeff,outlier,outs,opts] = am_solver(Xdot, phiX, opts)
%% Parameters and Defaults
switch nargin % if number of inputs is 2 then use default parameters
    case 2
        opts =[];
end

if isfield(opts,'maxit'), maxit = opts.maxit; else maxit = 100; end
if isfield(opts,'tol'), tol = opts.tol; else tol = 1e-2; end
if isfield(opts,'coef_thres'), coef_thres = opts.coef_thres; else coef_thres = 0.1; end
if isfield(opts,'mu'), mu = opts.mu; else mu = 0.0125; end
    
%% Operators    
% Shrink2 operators soln = argmin gamma/2 ||d-w||_2^2 + sum_i||d^i||_2,
% where d^i denotes the ith row of d
shrink2 = @(w,gamma) repmat(max(1-1./(gamma*sqrt(sum(w.^2,2))),0), 1, size(w,2)).*w;

%% Algorithm
% Initialization
outlier = zeros(size(Xdot));
b = randn(size(Xdot));
coeff = phiX\(Xdot - outlier - b );

iter = 0;
errord = 1.0; % initialize the difference of two consecutive values of outlier iterates
while ((iter <= maxit) && (errord > tol) )
    iter = iter+1;
    coeffs_old = coeff;
   
    % C-subproblem: compute the coefficients
    coeff = phiX\(Xdot - outlier - b );
    for iter_inner = 1:5 % do hard-thresholding to ensure the sparsity of C
        smallinds = (abs(coeff)<coef_thres); % find small coefficients
        coeff(smallinds) = 0; % hard-thresholding
        for ind = 1:size(Xdot,2)
            biginds = ~smallinds(:,ind);
            coeff(biginds,ind) = phiX(:,biginds)\(Xdot(:,ind) - outlier(:,ind) - b(:,ind));
        end
    end
    relerror_coeff(iter) = norm(coeffs_old - coeff,'fro');
    
    % E-subproblem: compute the outliers
    Eold = outlier;
    outlier = shrink2(Xdot - phiX*coeff - b , mu);
    
    errord = norm(outlier - Eold,'fro');
    
    relerror_outlier(iter) = norm(outlier-Eold,'fro');
    len_outlier(iter) = length(find(outlier(:,1)));
    % update multiplier
    b = b+ phiX*coeff  + outlier - Xdot;
end

%% Assign outputs
outs.iter = iter;
outs.relerror_coeff = relerror_coeff;
outs.relerror_outlier = relerror_outlier;
outs.len_outlier = len_outlier;

opts.maxit = maxit;
opts.tol = tol;
opts.coef_thres = coef_thres;
opts.mu = mu;

