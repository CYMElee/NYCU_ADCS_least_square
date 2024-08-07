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
record_theta3_hat = zeros(length(t),3);
record_theta = zeros(length(t),3);
record_Omega = zeros(length(t),3);
record_icl = zeros(length(t),3);

% ICL record 
record_y = zeros(3,3,length(t));
record_m_p = zeros(3,length(t));

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
   


    record_theta_hat(i,:) = Theta_hat';
    record_theta(i,:) = Theta';
    

    record_m_p(:,i) = M_p*dt; 
    record_y(:,:,i) = Y_cl*dt;
    

    M_p = -kr*eR - kw*eW - Y_j*Theta_hat;

    M =[M_p;0];

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
   

end

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
 [g_hat_map] = hat_map();
 [omega_hat_map] = hat_map();

 %/*************************************************/
 %/                     Psi matrix                  /
 %/*************************************************/

 % Psi matrix is a 3N*9 matrix %
 %(1):Integral the centripetal force %

 %(2):Integral gravitational acceleration(skew-symmetric form) %


 %(3):using (1)(2),obtain Psi matrix %
 Psi =[omega_ab-omega_ab_init+ , ];
 Psi_N =[Psi_N;Psi];



 Z_N = 



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









