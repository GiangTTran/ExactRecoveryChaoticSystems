%==========================================================================
% dictionary4: dictionary matrix built from the 4D data
%          
% Input: 
%       U(m x 4): time-varying measurements [x y z w], kth row is the measurement value at time k*dt 
%                         where 4 = dimension of the ODE system
%                               m = # measurements
%                               r = # basis functions
% Output:
%       phiX(m x r): dictionary matrix [1 x y z w x^2 xy xz y^2 ... w^p]
%                         where x,y,z,w... are column vectors
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
function phiX = dictionary4(U)

% phiX =[1 X Y Z X^2 XY XZ Y^2 YZ Z^2 ... ]
phiX = zeros(size(U,1),70);
% 1 - 1 column
phiX(:,1) = ones(size(U,1),1);% 1
% X - 4 columns
phiX(:,2:5) = U;               % x y z w
% X^2 - 10 columns
phiX(:,6:9) = repmat(U(:,1),1,4).*U(:,1:4); % x^2, xy, xz, xw
phiX(:,10:12) = repmat(U(:,2),1,3).*U(:,2:4);  % y^2, yz, yw
phiX(:,13:14) = repmat(U(:,3),1,2).*U(:,3:4); % z^2, zw
phiX(:,15) = U(:,4).^2; % w^2
% X^3 - 20 columns
phiX(:,16:25) = repmat(U(:,1),1,10).*phiX(:,6:15); % x^3, x^2y, x^2z, x^2w,xy^2, xyz, xyw, xz^2,xzw, xw^2
phiX(:,26:31) = repmat(U(:,2),1,6).*phiX(:,10:15); %y^3, y^2z, y^2w, yz^2, yzw, yw^2
phiX(:,32:34) = repmat(U(:,3),1,3).*phiX(:,13:15); % z^3, z^2w, zw^2
phiX(:,35) = U(:,4).^3; % w^3
% % % X^4 - 35 columns
phiX(:,36:55) = repmat(U(:,1),1,20 ).*phiX(:,16:35); % x^4,...
phiX(:,56:65) = repmat(U(:,2),1,10).*phiX(:,26:35); % y^4, ...
phiX(:,66:69) = repmat(U(:,3),1,4).*phiX(:,32:35); % z^4, z^3w, z^2w^2, zw^3
phiX(:,70) = U(:,4).^4; % w^4
return
