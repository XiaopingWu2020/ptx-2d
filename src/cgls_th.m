function [x, nIters,rho,eta] = cgls_th (A, b, th)
% CGLS_TH Solve least square problems based on a threshold
% modified from cgls.m by Hansen. The iteration stops if the difference between 2
% successive residual norms is less than a pre-defined threshold.
% 
% CGLS Conjugate gradient algorithm applied implicitly to the normal equations.
%
% Performs k steps of the conjugate gradient algorithm applied
% implicitly to the normal equations A'*A*x = A'*b.
%
% The routine returns all k solutions, stored as columns of
% the matrix X.  The corresponding solution and residual norms
% are returned in the vectors eta and rho, respectively.
%
% If the singular values s are also provided, cgls computes the
% filter factors associated with each step and stores them
% columnwise in the matrix F.
%
% Reorthogonalization of the normal equation residual vectors
% A'*(A*X(:,i)-b) is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.
% 
% References: A. Bjorck, "Numerical Methods for Least Squares Problems",
% SIAM, Philadelphia, 1996.
% C. R. Vogel, "Solving ill-conditioned linear systems using the
% conjugate gradient method", Report, Dept. of Mathematical
% Sciences, Montana State University, 1987.
% 
% Per Christian Hansen, IMM, 07/02/97.
%
% Usage: [x, nIters,rho,eta] = cgls_th (A, b, th)
%
% Returns
% -------
% x: solution
% nIters: num of iterations
%
% rho: residual err, |Ax - b|
%
% eta: norm(x)
%
%
% Expects
% -------
% A: system matrix
% b: target vector
% th: threshold used to stop iterations. 
%
%
% See also: cgls
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Fri Jun 13 09:33:58 2008
%

% Initialization.
[~,n] = size(A); 

% Prepare for CG iteration.
x = zeros(n,1);
d = A'*b;
r = b;
rho = norm(r);
normr2 = d'*d;

% Update x and r vectors once.
Ad = A*d; alpha = normr2/(Ad'*Ad);
x  = x + alpha*d;
r  = r - alpha*Ad;
s  = A'*r;
rho_new = norm(r);
% Update d vector once.
normr2_new = s'*s;
beta = normr2_new/normr2;
normr2 = normr2_new;
d = s + beta*d;

nIters = 1;
% Iterate.
while (rho - rho_new) > th
  nIters = nIters + 1;
  
  rho = rho_new;
% Update x and r vectors.
  Ad = A*d; alpha = normr2/(Ad'*Ad);
  x  = x + alpha*d;
  r  = r - alpha*Ad;
  s  = A'*r;
  rho_new = norm(r);
% Update d vector.
  normr2_new = s'*s;
  beta = normr2_new/normr2;
  normr2 = normr2_new;
  d = s + beta*d;
end

% Compute norms, if required.
  if (nargout>2), rho = rho_new; end
  if (nargout>3), eta = norm(x); end
  
disp('-> Done...')
