

%% ADCS Reaction wheels control system simulation %%
clc; clear;
%% Plant parameters %
% The wheel's moment of inertia include the roter and wheel.
J_RW_testbed = 2.27373e-4; %(kg*m^2)
% All of the body frame's moment of inertia
J_AB_testbed = [43.249 0.6530 0.0340;...
                0 78.66 0.2628;...
                0 0 98.06;];



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
r_CM = [0.1722;0.0089;-0.1516];
J = [43.249,78.66,98.06,0.6530,0.0340,0.2628]';
m_sys = 24; %units:kg
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







