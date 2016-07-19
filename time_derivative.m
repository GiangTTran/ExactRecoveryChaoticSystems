%==========================================================================
% time_derivative: 1st/2nd numerical approximation of time derivative 
%       1st-order approximation:
%                               roc(x(t)) = (x(t+dt) - x(t))/dt
%       2nd-order approximation:
%                               roc(x(t)) = (x(t+dt) - x(t-dt))/(2*dt)
%
% Input: 
%       X(mxn): recorded data, kth row is the measurement value at time k*dt 
%       dt:     time step
%                    where m = number of measurements
%                          n = dimension of the ODE system
% Output: 
%       roc((m-1)xn: 1st-order approximation of time derivative
%       roc((m-2)xn: 2nd-order approximation of time derivative
%
% Note: The code can be generalized for non-equal time distribution. For
% example, the 2nd-order approximation will be 
%       roc(x(t_k)) = (x(t_{k+1}) - x(t_{k-1})) / (t_{k+1} - t_{k-1});
%
% Authors: Giang Tran and Rachel Ward
% Institution: The University of Texas at Austin
% Version: 1.0, Year: 2016
% More information can be found at: 
%     G. Tran and R. Ward, "Exact Recovery of Chaotic Sysmtems from Highly
%     Corrupted Data", https://arxiv.org/abs/1607.01067 
%==========================================================================
function roc = time_derivative(X,dt,type)
[m,n] = size(X);
if (type==1)
    roc = zeros(m-1,n);
    roc = (X(2:m,:) - X(1:m-1,:))/ dt; % forward Euler
else if (type==2)
        roc = zeros(m-2,n);
        roc(1:m-2,:) = (X(3:m,1:n) - X(1:m-2,1:n)) /(2*dt); % note roc(k,:) is the time derivative at time t_{k+1} 
    end
end
return