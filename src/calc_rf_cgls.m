function [rf,b] = calc_rf_cgls (sysmat, m, nchs, timewin, th)

% CALC_RF_CGLS Calculate rf pulses with conjugate gradient method.
%
% Usage: rf = calc_rf_cgls (sysmat, m, nchs, timewin, th)
%
% Returns
% -------
% rf: a nchs x ntimep matrix containing rf pulses
%
% Expects
% -------
% sysmat: system matrix A
% m: target vector
% nchs: num of Tx channels
% timewin: time window. defaults to [] (empty)
% th: threshold determining the termination of iterations. defaults to 5e-3
%
%
% See also: calc_rf_tikh calc_rf_addang cgls_th construct_target_vector
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Tue May 27 11:50:10 2008
%

if nargin < 5
  th = 5e-3;
end

if nargin < 4
  timewin = [];
end

[b,niters] = cgls_th(sysmat,m,th);
fprintf('-> Number of CG iterations = %d, |Ax - b| / |b| = %d\n',niters, ...
        norm(sysmat*b - m) / norm(m));

rf = reshape_rf(b,nchs,timewin);
