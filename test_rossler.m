%==========================================================================
% test the Rossler data with our optimization model
%      x' = -y-z; y' = x + alphaR*y ; z' = betaR + xz-gammaR*z
%
% Inputs:
%   parameters of the ODE 
%   params.: store all parameters 
%        dt: time step
%        Tfinal: final time
%        roc: (1 or 2) order to approximate time derivative
%
%      % (used for corrupted_data.m)
%        ratio_corrupted: corruption ratio 
%        sigma_corrupted: standard deviation of the Gaussian noise, added at the outlier(corrupted) location
%        min_blocklength: min length of the corrupted block
%        max_blocklength: max length of the corrupted block
%        sigma_noise: noise level added at every measurement, in addition to the outlier
%
%      % (used for am_solver.m)
%        tol   : tolerance for relative change of the outliers, default = 5e-3 
%        maxit : max number of iterations, default = 100
%        coef_thres: thresholding parameter for the coefficient (coeffs - hard thresholding)
%        mu    : weight of the term ||phiX*C + E - Xdot + b||_F^2
%
% Outputs: 
%    C: recovered coefficients
%    E: computed outliers
%    outs.: A struct with convergence information
%        iter  : #iterations
%        relerror_coeff: norm2 of two consecutive coeff's iterates each iteration
%        relerror_outlier: norm2 of two consecutive outlier's iterates each iteration
%        len_outlier: #nonzero entries in one column of outlier E each iteration
%
% Visualizations:
%    1. Plot Rossler system with outliers
%    2. Plot the true and recovered coefficients
%    3. Plot norm2 of two consecutive coeff's iterates vs Iterations
%    4. Plot norm2 of two consecutive outliers' iterates vs Iterations
%    5. Plot #nonzero entries in 1st column of outlier iterate vs Iterations
%    6. Print out the percentage of corruption (note: ratio_corrupted says about #corrupted block, each block has different length)
%    7. Print out the coefficient error (see definition in our paper) between the true and the computed coefficients
%    8. Print out whether the algorithm detects exactly the outliers
%    9. Print out number of iterations
%
% Authors: Giang Tran and Rachel Ward
% Institution: The University of Texas at Austin
% Version: 1.0, Year: 2016
% More information can be found at: 
%     G. Tran and R. Ward, "Exact Recovery of Chaotic Sysmtems from Highly
%     Corrupted Data", https://arxiv.org/abs/1607.01067 
%==========================================================================
close all; clear all; clc

%% Parameters
% solve the ODE using 4th order RK
params1.dt = 0.0005;
params1.Tfinal = 50.0;
params1.roc = 2; % 1st or 2nd order of time derivative approximation 

% parameters for generating corrupted data
params_data.ratio_corrupted = 0.004; % corruption ratio
params_data.sigma_corrupted = 50*params1.dt; % standard deviation of the Gaussian noise, added at the outlier(corrupted) location
params_data.min_blocklength = 5; % min length of the corrupted block
params_data.max_blocklength = 50; % max length of the corrupted block
params_data.sigma_noise = 0; % noise level added at every measurement, in addition to the outlier

% parameters for algorithm
params_alg.tol = 1e-5; % tolerance for relative change of the outliers, default = 5e-3
params_alg.maxit = 500; % max number of iterations, default = 100
params_alg.coef_thres = 0.05; % thresholding parameter for the coefficient (coeffs - hard thresholding)
params_alg.mu = 0.025; % weight of the term ||phiX*C + E - Xdot + b||_F^2

% Parameters of the Rossler system
alphaR = 0.2;
betaR = 0.2;% 470/19.0;%28
gammaR = 5.7;
U0 = [0.0 -6.78 0.02]; % Starting point from Rossler's paper - An Equation for Continuous Chaos

% Generate clean data
myfunc = @(t,x) [-x(2)-x(3); x(1) + alphaR*x(2) ; betaR + x(1).*x(3)-gammaR*x(3)];
[T U_clean] = ode45(myfunc,[0:params1.dt:params1.Tfinal],U0,odeset('RelTol',1e-12,'AbsTol',[1e-12,1e-12,1e-12]));

%% Generate Corrupted Data, Time Derivative of Data, and Dictionary 

[U,index_mislead,opts_data] = corrupted_data(U_clean,params1.dt,params_data);

% Approximate time derivative of data
Xdot = time_derivative(U,params1.dt,params1.roc); 

% Build Dictionary
phiX = dictionary3(U);

% matching the dimension of U,Xdot and phiX
U(end,:) = []; % associate with forward Euler of time derivative
phiX(end,:) = []; % associate with forward Euler of time derivative

if (params1.roc == 2)
    U(1,:)=[];
    phiX(1,:)= [];
end

%% Main Algorithm
%   C: recovered coefficients
%   E: numerical outliers
%   outs.iter: #iterations
%   outs.error_coeff: norm2 of two consecutive coeff's iterates each iteration
%   outs.error_outlier: norm2 of two consecutive outlier's iterates each iteration
%   outs.len_outlier: #nonzero entries in one column of outlier E each iteration

[C,E,outs,opts_alg] = am_solver(Xdot, phiX,params_alg);

names = [fieldnames(params1);fieldnames(opts_data); fieldnames(opts_alg)];
params = cell2struct([struct2cell(params1);struct2cell(opts_data); struct2cell(opts_alg)],names, 1)
%% Visualization
% Visualize the sample data
figure('name','Rossler with corrupted data');
plot3(U(:,1), U(:,2), U(:,3),'.','Linewidth',1.25); 
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('z','FontSize',18)
set(gca,'FontSize',18)
axis tight; grid on; view(-40,55)


% Compare the recovered vs the true coefficients
% True Solution for the coefficient C
%% True Solution
a0 = [0, 0, -1, -1,0,0,0]';%[1, x, y, z, x^2, xy, xz]
b0 = [0, 1, alphaR, 0,0,0,0]';
c0 = [betaR,0,0,-gammaR, 0, 0,1]';

coeff_true = [a0 b0 c0];
C_true= [coeff_true; zeros(size(phiX,2)-7,3)];

figure('name','True and Recovered Coefficients')
    subplot(311); plot(C(:,1),'r*'); hold on; plot(C_true(:,1),'+b');
    legend('Numerical Soln','Ground-truth Soln')
    subplot(312); plot(C(:,2),'r*');hold on; plot(C_true(:,2),'+b');
    subplot(313); plot(C(:,3),'r*');hold on; plot(C_true(:,3),'+b');
    
% Plot errors vs Iterations
figure('name','Errors vs Iterations'); 
    subplot(311);semilogy(outs.relerror_coeff,'*'); title('Relative Error of Coefficients vs Iterations')
    subplot(312); plot(outs.len_outlier,'*'); title('Length of Outlier Variable vs Iterations')
    subplot(313); semilogy(outs.relerror_outlier,'*'); title('Error of Outlier Variable vs Iterations')

% Percentage of corruption
fprintf('The percentage of block corrupting data %2.2f\n',100*length(index_mislead)/size(U_clean,1))
fprintf('Number of Iterations is %3d\n',outs.iter)
% Compute the coefficient error (see definition in our paper)

error_coeffMax = (C(abs(C)>0) - C_true(abs(C)>0))./C_true(abs(C)>0);
fprintf('The percentage of error coefficient %2.4f\n',max(abs(error_coeffMax))*100)

%% Locations of outliers
% If we use the 2nd-order approximation:
    % note the system start from t_2->t_{m-1}. 
    % If the misled indices are from time t_48 to t_100, the mislead in derivative
    % is from t_47 to t_101. In the new system, the misled indices are from 46 to 100
% If we use the forward Euler - 1st order approximation
    % note the system start from t_1->t_{m-1}. 
    % If the misled indices are from time t_48 to t_100, the mislead in derivative
    % is from t_47 to t_100. In the new system, the misled indices are from 47 to 100
if (params1.roc==2)
    tmp2 = index_mislead-2;
else 
    tmp2 = index_mislead-1;
end

index_true =unique([index_mislead tmp2]);

index_computed = find(E(:,1)); % expect index_computed = index_true
if (isempty(setdiff(index_computed,index_true)) && isempty(setdiff(index_true,index_computed)))
    fprintf('The algorithm detects exactly the locations of outliers \n')
else if ((length(setdiff(index_computed,index_true)) + length(setdiff(index_true,index_computed))) <=2)
        fprintf('The algorithm detects the locations of outliers except at most 2 locations \n')
    else fprintf('The algorithm doesnot detect well enough the locations of outliers\n')
    end
end

    
    




