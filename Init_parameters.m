

%% ADCS Reaction wheels control system simulation %%
clc; clear;
%% Plant parameters %
% The wheel's moment of inertia include the roter and wheel.
J_RW_testbed = 2.27373e-4; %(kg*m^2)
% All of the body frame's moment of inertia
J_AB_testbed = [5.47 0 0;...
                0 5.47 0;...
                0 0 5.47];



%gravitational acceleration
g = 9.81; %kg/m^2
g_body = [0,0,0]'; %3x1 matrix
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

ext_Torque = [0;0;0];
Rq =[1;0;0;0];
q0 =[1;0;0;0];

omega_ab_init = [0,0,0]'; % using for Psi matrix
omega_ab = [0,0,0]';
omega_ab_prev = [0,0,0]';
omega_dot_ab = [0,0,0]';


omega_mo = [0,0,0,0]';
omega_dot_mo = [0,0,0,0]';
omega_mo_prev = [0,0,0,0]';



r = [0,0,0]; 

%control input
M_p = [0,0,0]';   
M = [0,0,0]';

Z_N = [0,0,0]';
GF =[0,0,0];
CF =zeros(3,6); 
CI = [0,0,0]';
% Real values
r_CM = [0.001;0.002;0.003];
r_CMXM = [60.56*0.001,60.56*0.002,60.56*0.003];
J_MOI = [5.47,5.47,5.47]';
J_POI = [0,0,0]';
m_sys = 60.56; %units:kg
%Estimate values
r_CM_hat  = [0;0;0];
J_hat = [0;0;0;0;0;0];

% Robust control control gains Gp Gr Gamma
Gp = [25.2,0,0; ...
      0,25.2,0; ...
      0,0,28.8];

Gr = [20,0,0; ...
      0,20,0; ...
      0,0,32];

Gamma = 20;







