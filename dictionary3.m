%==========================================================================
% dictionary3: dictionary matrix built from the 3D data
%          
% Input: 
%       U(m x 3): time-varying measurements [x y z], kth row is the measurement value at time k*dt 
%                         where 3 = dimension of the ODE system
%                               m = # measurements
%                               r = # basis functions
% Output:
%       phiX(m x r): dictionary matrix [1 x y z x^2 xy xz y^2 ... z^p]
%                         where x,y,z,... are column vectors
%
% Note: max degree of the multivariate polynomial is 4
%
% Authors: Giang Tran and Rachel Ward
% Institution: The University of Texas at Austin
% Version: 1.0, Year: 2016
% More information can be found at: 
%     G. Tran and R. Ward, "Exact Recovery of Chaotic Sysmtems from Highly
%     Corrupted Data", https://arxiv.org/abs/1607.01067 
%==========================================================================
function phiX = dictionary3(U)

phiX = zeros(size(U,1),35);
% 1 - 1 column
phiX(:,1) = ones(size(U,1),1);% 1
% X - 3 columns
phiX(:,2:4) = U;               % x y z
% X^2 - 6 columns
phiX(:,5:7) = repmat(U(:,1),1,3).*U(:,1:3); % x^2, xy, xz
phiX(:,8:9) = repmat(U(:,2),1,2).*U(:,2:3);  % y^2, yz
phiX(:,10) = U(:,3).^2; % z^2
% X^3 - 10 columns
phiX(:,11:16) = repmat(U(:,1),1,6).*phiX(:,5:10); % x^3, x^2y, x^2z, xy^2, xyz, xz^2
phiX(:,17:19) = repmat(U(:,2),1,3).*phiX(:,8:10); %y^3, y^2z, yz^2
phiX(:,20) = U(:,3).^3; % z^3
% X^4 - 15 columns
phiX(:,21:30) = repmat(U(:,1),1,10 ).*phiX(:,11:20); % x^4, x^3y, x^3z, x^2y^2, x^2yz, x^2z^2, xy^3, xy^2z, xyz^2, xz^3
phiX(:,31:34) = repmat(U(:,2),1,4).*phiX(:,17:20); % y^4, y^3z, y^2z^2, yz^3
phiX(:,35) = U(:,3).^4; % z^4
% % X^5 - 21 columns
% phiX(:,36:50) = repmat(U(:,1),1,15).*phiX(:,21:35); % U(:,1).*X^4 - 15 columns
% phiX(:,51:55) = repmat(U(:,2),1,5).*phiX(:,31:35); % y^5, y^4z, y^3z^2, y^2z^3, yz^4 - 5 columns
% phiX(:,56) = U(:,3).^5; % z^5
return
