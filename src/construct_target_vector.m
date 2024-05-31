function targetv = construct_target_vector (targetMxyPattern, spatialMask, ...
                                                  fa, ttype)

% CONSTRUCT_TARGET_VECTOR construct target vector for Transmit SENSE pulse
% design.
%
% Usage: targetv = construct_target_vector (targetMxyPattern, spatialMask, fa, ttype)
%
% Returns
% -------
% targetv: 
%
% Expects
% -------
% 
% targetMxyPattern: should be complex valued pattern with nonzeros for
% excitation (a.u.). when empty, the spatialMask will be taken as
% targetMxyPattern.
% 
% spatialMask: logic
% 
% fa: max nominal flip angle in deg. defaults to 5 deg.
% 
% ttype: a char indicating what kind of target is specified. can be 's' for
% small tip angle pulse design where target is defined by Mxy/M0, or 'l' for
% large tip angle where target is defined by flip angles, which can also be used
% as complex flip angle target (refer to Boulant and Hoult MRM 2012 for more
% details on complex flip angles). defaults to 's'.
% 
% 

  if isempty(targetMxyPattern)
    targetMxyPattern = spatialMask;
  end

  if ~isnumeric(targetMxyPattern)
    targetMxyPattern = double(targetMxyPattern);
  end
  
  if isnumeric(spatialMask)
    spatialMask = logical(spatialMask);
  end
  
  if nargin < 3
    fa = 5; % nominal flip angle in deg.
  end 
  
  if nargin < 4
    ttype = 's';
  end

% normalization
  targetMxyPattern = targetMxyPattern ./ max(targetMxyPattern(:));
      
  switch ttype
   case 's'
    targetv = sin(fa/180*pi)* targetMxyPattern(spatialMask);
   case 'l'
    targetv = (fa/180*pi)* targetMxyPattern(spatialMask);
   otherwise
  end  
  
  disp('=> Target vector constructed.');
