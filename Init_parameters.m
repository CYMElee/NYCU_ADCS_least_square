

%% ADCS Reaction wheels control system simulation %%
clc; clear;
%% Plant parameters %
% The wheel's moment of inertia include the roter and wheel.
J_RW_testbed = 2.27373e-4; %(kg*m^2)
% All of the body frame's moment of inertia
J_AB_testbed = [6.2 0 0;...
                0 6.2 0;...
                0 0 5.7;];

%gravitational acceleration
g = 9.81; %kg/m^2
g_inertial = [0,0,g]';
g_body = [0,0,0]';
% The direction consine matrix(DCM) from RWs to Body frame
A_w = [1/2 -1/2 -1/2 1/2;
        1/2 1/2 -1/2 -1/2;
        1/(2)^0.5 1/(2)^0.5 1/(2)^0.5 1/(2)^0.5];


% H_w matrix is the A_w after power optimization.
H_w = [1/2 -1/2 -1/2 1/2;
       1/2 1/2 -1/2 -1/2;
       1/(2)^0.5 1/(2)^0.5 1/(2)^0.5 1/(2)^0.5;
       1 -1 1 -1];

% Init Condition %

ext_Torque = [0,0,0];
Rq =[1;0;0;0];
q0 =[1;0;0;0];
%Initialization for function "Platform_Dynamics_Func"%
omega_ab = [0,0,0]';
omega_ab_prev = [0,0,0]';
omega_dot_ab = [0,0,0]';
omega_dot_ab_prev =[0,0,0]';

omega_mo = [0,0,0,0]';
omega_dot_mo = [0,0,0,0]';
omega_mo_prev = [0,0,0,0]';
omega_dot_mo_prev = [0,0,0,0]';

M_prev = [0,0,0,0]';
r = [0,0,0]; 

%control input
M_p = [0,0,0]';   
M = [0,0,0]';


% Real values
r_CM = [1,2,3];
J = [6.2,6.2,5.7]'; 
%Estimate values
r_CM_hat  = [0;0;0];
J_hat = [0;0;0];








