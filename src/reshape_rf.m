function rf = reshape_rf (b, nchs, timewin)
% RESHAPE_RF Reshape the solution for rf pulses. the solution comes from a certain optim problem.
%
% Usage: rf = reshape_rf (b, nchs, timewin)
%
% Returns
% -------
% rf: a nchs x ntimep matrix containing rf pulses
%
% Expects
% -------
% b: solution from a certain optim problem
% nchs: num of Tx channels
% timewin: time window. defaults to [] (empty)
%
%
% See also: 
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Tue May 27 10:34:53 2008
%

if nargin < 3
  timewin = [];
end

rf1 = reshape(b,[],nchs);

if isempty(timewin)
  rf = rf1.';
else
  rf = complex(zeros(nchs, length(timewin)));
  rf(:,logical(timewin)) = rf1.';
end
