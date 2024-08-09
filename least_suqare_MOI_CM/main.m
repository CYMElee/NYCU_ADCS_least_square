% "ADCS Testbed" Simlation Main file %
% Initialization 
clear;clc;

addpath("../");
Init_parameters;
% Simulation time setting
t = [1:1:30000];
dt = 0.001;

% array for recoding the state (for plot)
record_R = zeros(length(t),3);

record_Omega = zeros(length(t),3);


% Psi record
record_centripetal_force = zeros(3,3,length(t));

record_gravitational_acceleration = zeros(3,3,length(t));

% M record

record_M = zeros(length(t),3);


for i=1:length(t)

    %% Trans Gravitational acceleration from inertial frame to body fix frame
    if i==1
        g_body = g*[-sin(0);sin(0)*cos(0);cos(0)*cos(0)];
    else
        g_body = g*[-sin(record_R(i-1,2));sin(record_R(i-1,1))*cos(record_R(i-1,2));cos(record_R(i-1,1))*cos(record_R(i-1,2))];
    end
        ext_Torque = cross(r_CM,m_sys*g_body);
   
   

    %% control input
    
    % if Euler angle's x,y belong(-5,555) we assign torque directly
    a1 = sin(r(3)/2)*cos(r(2)/2)*cos(r(1)/2)-cos(r(3)/2)*sin(r(2)/2)*sin(r(1)/2);
    a2 = cos(r(3)/2)*sin(r(2)/2)*cso(r(1)/2)+sin(r(3)/2)*cos(r(2)/2)*sin(r(1)/2);
    a3 = cos(r(3)/2)*cos(r(2)/2)*sin(r(1)/2)-sin(r(3)/2)*sin*(r(2)/2)*cos(r(1)/2);
    alpha = [a1; ...
             a2; ...   
             a3];
    alpha_4 = cos(r(3)/2)*cos(r(2)/2)*cos(r(1)/2)+sin(r(3)/2)*sin(r(2)/2)*sin(r(1)/2);
    [alpha_hat] = hat_map(alpha);

    %else
    if  abs(r(3)) < 0.0872 && abs(r(2)) < 0.0872
        M_p = 2*[sin((2*pi*i*dt)/5);sin((2*pi*i*dt)/5);sin((2*pi*i*dt)/5)];
    else
        M_p = -0.5*(((alpha_hat)+alpha_4*eye(3))*Gp + Gamma*(1-alpha_4)*eye(3))*alpha-Gr*omega_ab;
    end
    record_M = M_p';
    M = [M_p;0];
    %% using A.B desire M to get R.W generate M
    omega_dot_mo = -(inv(H_w)/J_RW_testbed)*M;
    omega_mo = omega_mo_prev + omega_dot_mo*dt;
    omega_mo_prev = omega_mo; 

    %% system dynamics

    omega_dot_ab = inv(J_AB_testbed)*...
                       (ext_Torque-((A_w*J_RW_testbed*omega_dot_mo)+...
                        cross(omega_ab_prev, (A_w*J_RW_testbed*omega_mo+J_AB_testbed*omega_ab_prev))));

    omega_ab = omega_ab_prev + omega_dot_ab*dt;
    omega_ab_prev = omega_ab;
   
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
 %/  Define angular rate matrix( 3x6 matrix) /%
 %/******************************************/%
 [omega_N] = Angular_rate_matrix(omega_ab); % angular rate matrix for OMEGA(N)
 [omega_init] = Angular_rate_matrix(omega_ab_init); % angular rate matrix for OMEGA(t0)
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
 Psi_temp = [omega_N-omega_init+CF ,GF ];
 Psi_N = [Psi_N;Psi_temp];

 % The pseudo-inverse of Psi_N

 inv_Psi_N = inv(Psi_N'*Psi_N)*Psi_N';

 %/*************************/%
 %/            Z_N          /%
 %/*************************/%

% (1): control input integral
if i == 1
    Z_N = record_M(N,:)'*dt;
else
    for N= 1:length(i)
        CI =+ (record_M(N,:)'*dt); %CI is 1x3 matrix    
    end
   Z_N = [Z_N;CI]; %Z_N is 3Nx1 matrix, in here,CI should be transport.
end
 


% using least square to estimate X
 X_hat = Psi_N'*inv_Psi_N*Psi_N'*Z_N;

 J_hat = X_hat(1:6,:);
 r_CM_hat = X_hat(7:9,:);

 end

%% Plot


%% Omega plot

figure;
ax1 = nexttile;
plot(ax1, ...  
          t, record_Omega(1:length(record_Omega),1),'--'...
          );
title("Angular rate(rad/s)", FontSize=16);
xlabel("Time(ms)", FontSize=13);
ylabel("omega(rad/s)", FontSize=13);
legend("x");












