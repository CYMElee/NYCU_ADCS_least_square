  % "ADCS Testbed" Simlation Main file %
% Initialization 
clear;clc;

addpath("../");
Init_parameters;
% Simulation time setting
t =  [1:1:30000];
dt = 0.001;

% array for recoding the state (for plot)
record_R = zeros(length(t),3);
record_eR =zeros(length(t),3);
record_Omega = zeros(length(t),3);


% Psi record
record_centripetal_force = zeros(3,3,length(t));

record_gravitational_acceleration = zeros(3,3,length(t));

% M record

record_M = zeros(length(t),3);

%desire Attitude and desire Omega(platform bodyframe)
r = [0,0,0];



for i=1:length(t)
 
    % set the desire Attitude,Angular velocity, Angular acceleration
    rd = [0,0,sin(i*dt)];%xyz 
    Od_raw = [0,0,cos(i*dt)];
    Odd_raw = [0,0,-sin(i*dt)];


     % let the trajectory compare with the Euler angle principle
     %Reference:https://ocw.mit.edu/courses/2-154-maneuvering-and-control-of-surface-and-underwater-vehicles-13-49-fall-2004/bc67e15b31b4f30aceabef2a66a6229d_lec1.pdf
     [Od,Odd]=Get_Omega_Desire(rd,Od_raw,Odd_raw,i,dt);
     Od = Od';
     Odd = Odd';
   

    record_Rd(i,:) = rd;
    record_Od(i,:) = Od;
    record_theta(i,:) = Theta.';
    record_theta_hat(i,:) = Theta_hat.';
    %% Get the tracking error for Attitude and Omega
    rd = [rd(3),rd(2),rd(1)];
    Rd = eul2rotm(rd);
    R = eul2rotm(r);
    [eR] = Attitude_Error(Rd,R);
    [eW] = Omega_Error(omega_ab_prev,Od,R,Rd);
    
    [o_b] = omega_bar(omega_ab_prev,R,Rd,Od,Odd);
    record_eR(i,:) = [eR(1),eR(2),eR(3)];
    record_eW(i,:) = [eW(1),eW(2),eW(3)];


    %% Trans Gravitational acceleration from inertial frame to body fix frame
    g_body = g_inertial*[-sin(record_R(i-1,2)),sin(record_R(i-1,1))*cos(record_R(i-1,2)),cos(record_R(i-1,1))*cos(record_R(i-1,2))];
    ext_Torque = cross(r_CM,m*g_body);
    ext_Torque = ext_Torque';
   

    %% using A.B desire M to get R.W generate M
    omega_dot_mo = -(inv(H_w)/J_RW_testbed)*M;
    omega_mo = omega_mo_prev + omega_dot_mo*dt;
    omega_dot_mo_prev = omega_dot_mo;
    omega_mo_prev = omega_mo;

    %% system dynamics

    omega_dot_ab = inv(J_AB_testbed)*...
                       (ext_Torque-((A_w*J_RW_testbed*omega_dot_mo)+...
                        cross(omega_ab_prev, (A_w*J_RW_testbed*omega_mo+J_AB_testbed*omega_ab_prev))));

    omega_ab = omega_ab_prev + omega_dot_ab*dt;
    omega_ab_prev = omega_ab;
    omega_dot_ab_prev = omega_dot_ab;

    %% using dynamic get Attitude(Quaternion)

    omega_x =omega_ab(1);
    omega_y =omega_ab(2);
    omega_z =omega_ab(3);
    omega_Matrix = [0, -omega_x, -omega_y, -omega_z;...
                    omega_x, 0, omega_z, -omega_y; ...
                    omega_y, -omega_z, 0, omega_x;...
                    omega_z, omega_y, -omega_x, 0];
    omega_Abs = norm(omega_ab);

    if omega_Abs == 0  %prevent devided by 0 occur when omega is 0
        Rq = q0;
    else
        I4 = eye(4);
        Quaternion = (cos((omega_Abs*dt)/2)*I4+((1/omega_Abs)*sin((omega_Abs*dt)/2)*omega_Matrix))*q0;
        q_norm = norm(Quaternion);
        Rq = (Quaternion/q_norm);
    end 
    q0 = Rq;

    % get Attitude(Euler angle)
    r = quat2eul(Rq.');
    %% using dynamic get Attitude(Euler angle)

    record_R(i,:) =[r(3),r(2),r(1)];
 
    record_Omega(i,:) =[omega_ab(1),omega_ab(2),omega_ab(3)];
   



 %% least square
 %/******************************************/%
 %/        Define angular rate matrix        /%
 %/******************************************/%
 Angular_rate_matrix = [omega_ab(1),0,0,omega_ab(2),omega_ab(3),0; ...
                        0,omega_ab(2),0,omega_ab(1),0,omega_ab(3); ...
                        0,0,omega_ab(3),0,omega_ab(1),omega_ab(3)];
 %/*******************************************/%
 %/           Skew-symmetric matrix           /%
 %/*******************************************/%
 [g_hat_map] = hat_map(g_body);
 [omega_hat_map] = hat_map(omega_ab);
 record_centripetal_force(:,:,i) =  omega_hat_map*Angular_rate_matrix;  
 record_gravitational_acceleration(:,:,i) = g_hat_map;
 %/*************************************************/
 %/                     Psi matrix                  /
 %/*************************************************/

 % Psi matrix is a 3N*9 matrix %
 %(1):Integral the centripetal force %

 for N= 1:length(i)
    CF =+ (record_centripetal_force(:,:,N))*dt;
 end

 %(2):Integral gravitational acceleration(skew-symmetric form) %
 for N= 1:length(i)
     GF =+ (record_gravitational_acceleration(:,:,N)*dt);    
 end
 %(3):using (1)(2),obtain Psi matrix %
 Psi_temp =[omega_ab-omega_ab_init+CF ,GF ];
 Psi_N = [Psi_N;Psi_temp];

 %/*************************/%
 %/            Z_N          /%
 %/*************************/%

% (1): control input integral
 for N= 1:length(i)
     CI =+ (record_M(N,:)*dt); %CI is 1x3 matrix    
 end
 Z_N = [Z_N;CI']; %Z_N is 3Nx1 matrix, in here,CI should be transport.

% using least square to estimate x

 X_hat = Psi_N'*inv(Psi_N)*Psi_N'*Z_N;

 
 end

%% Plot

%% Attitude Error 
figure;
ax1 = nexttile;
plot(ax1, t, record_eR(1:length(record_eR),3)...
          );
title("Attitude Error(z)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Euler Angle(rad)", FontSize=13);
legend("z");

%% Omega plot

%

figure;
ax2 = nexttile;
plot(ax2, ...  
          t, record_Od(1:length(record_Od),3),...
          t, record_Omega(1:length(record_Omega),3),'--'...
          );
title("Omgea", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Euler Angle(rad/s)", FontSize=13);
legend("zd","z");


%% Attitude plot


figure;
ax3 = nexttile;
plot(ax3,...
          t, record_Rd(1:length(record_Rd),3),...
          t, record_R(1:length(record_R),3),'--'...
          );
title("Attitude(Z)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Euler Angle(rad)", FontSize=13);
legend("zd","z");

figure;
ax4 = nexttile;
plot(ax4,...
          t, record_Rd(1:length(record_Rd),2),...
          t, record_R(1:length(record_R),2),'--'...
          );
title("Attitude(Y)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Euler Angle(rad)", FontSize=13);
legend("yd","y");



figure;
ax5 = nexttile;
plot(ax5,...
          t, record_Rd(1:length(record_Rd),1),...
          t, record_R(1:length(record_R),1),'--'...
          );
title("Attitude(X)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Euler Angle(rad)", FontSize=13);
legend("xd","x");



%% Estimate 
figure;
ax6 = nexttile;
plot(ax6, ...
          t, record_theta(1:length(record_theta),3),...
          t, record_theta_hat(1:length(record_theta_hat),3),"--"...
          );
title("Estimates of the Moment of Inertia(Izz)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Moment of Inertia(kg*m^2)", FontSize=13);
legend("Izz","Izz estimated");

figure;
ax7 = nexttile;
plot(ax7, ...
          t, record_theta(1:length(record_theta),2),...
          t, record_theta_hat(1:length(record_theta_hat),2),"--"...
          );
title("Estimates of the Moment of Inertia(Iyy)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Moment of Inertia(kg*m^2)", FontSize=13);
legend("Iyy","Iyy estimated");

figure;
ax8 = nexttile;
plot(ax8, ...
          t, record_theta(1:length(record_theta),1),...
          t, record_theta_hat(1:length(record_theta_hat),1),"--"...
          );
title("Estimates of the Moment of Inertia(Ixx)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("Moment of Inertia(kg*m^2)", FontSize=13);
legend("Ixx","Ixx estimated");









