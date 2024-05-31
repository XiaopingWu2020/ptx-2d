%%% pTx2D.m --- a script describing how to design 2D ptx pulses
% 
% Filename: pTx2D.m
% Description: 
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Maintainer: 

clear
addpath('./src/');
addpath('./bloch_sim/');

%% load calibration data etc.
load   data.mat 
% data.mat contains 16-ch B1+ maps, spatial mask, B0 map, 
% 2D gradient, gradient dwell time (in sec), field of excitation (fox, in m), 
% and an example target excitation pattern

%
sysmat = construct_system_matrix(grad,b1maps,mask,fox,[],[],[],dt);

targ = M;
nchs= size(b1maps,4);
rf = calc_rf_cgls(sysmat,construct_target_vector(targ,mask,10),nchs);

% bloch simulation
mxypatptx2d = run_bloch_sim (rf,grad,b1maps,mask,fox,b0map,0,[],dt);
figure, myimagesc(abs(mxypatptx2d),mask), colormap jet,axis equal


