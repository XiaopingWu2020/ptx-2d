function systemMatrix = construct_system_matrix (varargin)

% CONSTRUCT_SYSTEM_MATRIX Construct the system matrix for designing small tip
% angle pTx RF pulses in spatial domain.
%
% Usage: systemMatrix = construct_system_matrix (gradWaveform, b1plusMapArray, spatialMask,
% fieldOfExcitation, b0Map, temporalWindow, b0, dt, poffset, gamma)
%
% Returns
% -------
% systemMatrix: Ns-by-(Nt*Nc) matrix.
%
% Expects
% -------
% gradWaveform: 2,3-by-Nt array (tesla/m).
% b1plusMapArray: a 4D array [x,y,z,nCoils] with the last dimension for the number of coils.
% spatialMask: 2D or 3D logical array.
% fieldOfExcitation: 3-element vector [fox_x, fox_y, fox_z] (m).
% b0Map: 2D or 3D B0 inhomogeneity map (tesla) and its empty by default.
% temporalWindow: 1-by-Nt logical array indicating when to pulse and its empty by default.
% b0: a 1 x Nt vector for b0 eddy current, in tesla. defaults to [].
% dt: temporal resolution (sec). 4e-6 by default.
% 
% poffset: [offsetx,offsety,offsetz] in mm specifying the offset of FOV with
% respect to grad isocenter. this is mostly determined by RO Offset. defaults to
% [0 0 0].
% 
% gamma: gyromagnetic ratio (rad/sec/tesla). 2.675e8 by default for proton.
%
% 
% See also: construct_sysmat_lclta
% 

[grad, b1maps, mask, fox, b0map, timewin, b0, dt, poffset, method] = ...
    parseInputs(varargin{:});

  disp('=> Construction of the system matrix started...');
  
  kst = calc_ktraj_from_grad(grad,true,dt,method); % by gradient
  pha_main = calc_ktraj_from_grad(ones(1,length(timewin)),true,dt,method); % by B0 inhomogeneity
  pha_ec = calc_ktraj_from_grad(b0,true,dt,method); % by eddy current
  [b1arr, posarr] = create_array(b1maps,mask,fox,1e-3*poffset);
   
  phaarr = posarr*kst(:,timewin) + b0map(mask)*pha_main(:,timewin) + ...
           ones(size(b0map(mask)))*pha_ec(:,timewin);
       
  gamma = 2.675e8;
  systemMatrix = constrSysmat(phaarr,b1arr,dt,gamma);
  disp('=> Construction of the system matrix done.');
  
% ================
% 
% Local functions
% 
% ================

function [grad, b1maps, mask, fox, b0map, timewin, b0, dt,poffset,method] = ...
    parseInputs(varargin)

  error(nargchk(4, 10, nargin))

    % defaults.
  b0map = [];
  timewin = [];
  b0 = [];
  dt = 4e-6;
  poffset = [0, 0, 0];
  method = 's';

  switch nargin
   case 4
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
   case 5
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
      b0map = varargin{5};
   case 6
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
      b0map = varargin{5};
      timewin = varargin{6};
   case 7
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
      b0map = varargin{5};
      timewin = varargin{6};
      b0 = varargin{7};
   case 8
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
      b0map = varargin{5};
      timewin = varargin{6};
      b0 = varargin{7};
      dt = varargin{8};
   case 9
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
      b0map = varargin{5};
      timewin = varargin{6};
      b0 = varargin{7};
      dt = varargin{8};
      poffset = varargin{9};
   case 10
      grad = varargin{1};
      b1maps = varargin{2};
      mask = varargin{3};  
      fox = varargin{4};
      b0map = varargin{5};
      timewin = varargin{6};
      b0 = varargin{7};
      dt = varargin{8};
      poffset = varargin{9};
      method = varargin{10};
   otherwise
  end

  
  if size(grad,1) > 3
    grad = grad(1:3,:);
  end
  
  if ~islogical(mask)
    mask = logical(mask);
  end
  
  if isempty(b0map)
    b0map = zeros(size(mask));
  end

  if isempty(timewin)
    timewin = true(1,size(grad,2));
  end
  if ~islogical(timewin)
    timewin = logical(timewin);
  end

  if isempty(b0)
    b0 = zeros(1,size(grad,2));
  end

% check the validity
  if size(grad,1) ~= length(size(mask))
    error('The dimensions of grad and mask DONT match!')
  end
  if ndims(b1maps)~=4
    error('B1 map array should be 4D, i.e. [x,y,z,nCoils]!')
  end

  if size(b0map) ~= size(mask)
    error('The dimensions of b0map and mask DONT match!')
  end

  if length(timewin) ~= size(grad,2)
    error('The dimensions of timewin and grad DONT match!')
  end
% -----------------------

function systemMatrix = constrSysmat (phaArr, b1plusArray, dt, gamma)
% Returns
% -------
% systemMatrix: Ns-by-(Nt*Nc) matrix.
%
% Expects
% -------
% phaArr: Ns-by-Nt array for pseudo phases by gradient and B0 components (i.e.,
% main B0 inhomogeneity and b0 eddy current)
% b1plusArray: Ns-by-Nc array.
% dt: temporal resolution (sec). 4e-6 by default.
% gamma: gyromagnetic ratio (rad/sec/tesla). 2.675e8 by default for proton.

% construct the kernel matrix.
  disp('=> Construction of the kernel matrix started...');
  
  kernelMatrix = 1i .* gamma .* dt .* exp(1i .* phaArr);

  disp('=> Construction of the kernel matrix done.');
  % === the kernel matrix constructed ===

  % construct the system matrix.
  nCoils = size(b1plusArray,2);
  [ns,nt] = size(kernelMatrix);
  systemMatrix = complex(zeros(ns,nCoils*nt));
  for iCoil = 1 : nCoils,
    iIdx0 = (iCoil-1)*nt + 1;
    systemMatrix(:,iIdx0:(iIdx0+nt-1)) = diag(b1plusArray(:,iCoil)) * kernelMatrix;
  end
