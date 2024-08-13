% "ADCS Testbed" Simlation Main file %
% Initialization 
clear;clc;

addpath("../");
Init_parameters;
% Simulation time setting
t = [1:1:3000];
dt = 0.01;

% array for recoding the state (for plot)
record_R = zeros(length(t),3);
record_GF = zeros(3,3,length(t));

record_Omega = zeros(length(t),3);

record_MOI = zeros(length(t),3); % a length(t)*3 vector using store estimated MOI
record_POI = zeros(length(t),3); % a length(t)*3 vector using store estimated products of inertia

record_CM = zeros(length(t),3); % a length(t)*3 vector using store estimated C.M offset


% Psi record
record_centripetal_force = zeros(3,6,length(t));

record_gravitational_acceleration = zeros(3,3,length(t));

% M record

record_M = zeros(length(t),3);
record_Md = zeros(length(t),3);


for i=1:length(t)

    %% Trans Gravitational acceleration from inertial frame to body fix frame
    if i == 1
        g_body = g*[-sin(0);sin(0)*cos(0);cos(0)*cos(0)];
     
    else
        g_body = g*[-sin(record_R(i-1,2));sin(record_R(i-1,1))*cos(record_R(i-1,2));cos(record_R(i-1,1))*cos(record_R(i-1,2))];
    end
        ext_Torque = cross(r_CM,m_sys*g_body);
   
   

    %% control input
    
    % if Euler angle's x,y belong(-5,555) we assign torque directly
    a1 = sin(r(3)/2)*cos(r(2)/2)*cos(r(1)/2)-cos(r(3)/2)*sin(r(2)/2)*sin(r(1)/2);
    a2 = cos(r(3)/2)*sin(r(2)/2)*cos(r(1)/2)+sin(r(3)/2)*cos(r(2)/2)*sin(r(1)/2);
    a3 = cos(r(3)/2)*cos(r(2)/2)*sin(r(1)/2)-sin(r(3)/2)*sin(r(2)/2)*cos(r(1)/2);
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
    %record_M(i,:) = M_p';
    M = [M_p;0];
    %% using A.B desire M to get R.W generate M
    omega_dot_mo = -(inv(H_w)/J_RW_testbed)*M;
    omega_mo = omega_mo_prev + omega_dot_mo*dt;
    omega_mo_prev = omega_mo; 

    %% system dynamics
    % using R.C as momentum exchange devices.
    record_M(i,:) = - A_w*J_RW_testbed*omega_dot_mo -cross(omega_ab_prev,A_w*J_RW_testbed*omega_mo);

    omega_dot_ab = inv(J_AB_testbed)*...
                       (ext_Torque-(A_w*J_RW_testbed*omega_dot_mo)-...
                        cross(omega_ab_prev,A_w*J_RW_testbed*omega_mo)- ...
                        cross(omega_ab_prev, J_AB_testbed*omega_ab_prev));
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

    %% using quaternion integral to get Attitude

    if omega_Abs == 0  %prevent devided by 0 occur when omega is 0
       Rq = q0;
    else
        I4 = eye(4);
        Quaternion = (cos((omega_Abs*dt)/2)*I4+((1/omega_Abs)*sin((omega_Abs*dt)/2)*omega_Matrix))*q0;
        q_norm = norm(Quaternion);
        Rq = (Quaternion/q_norm);
    end 
    q0 = Rq;

   
    r = quat2eul(Rq.');


    


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
 record_centripetal_force(:,:,i) =  omega_hat_map*omega_N;  
 record_gravitational_acceleration(:,:,i) = g_hat_map;
 %/*************************************************/
 %/                     Psi matrix                  /
 %/*************************************************/

 % Psi matrix is a 3N*9 matrix %
 %(1):Integral the centripetal force %

 for N= 1:i
     CF_temp = record_centripetal_force(:,:,N)*dt;
     CF = CF + CF_temp;
 end

 %(2):Integral gravitational acceleration(skew-symmetric form) %
 for K= 1:i
     GF_temp = (record_gravitational_acceleration(:,:,K)*dt);
     GF = GF + GF_temp;
 end
 
 %(3):using (1)(2),obtain Psi matrix %
 if i == 1
    Psi_N = [omega_N-omega_init+CF ,GF ];
 else
     Psi_temp = [omega_N-omega_init+CF ,GF ];
     Psi_N = [Psi_N;Psi_temp];
 end
 record_GF(:,:,i) = GF;




 


 %/*************************/%
 %/            Z_N          /%
 %/*************************/%

% (1): control input integral
if i == 1
    Z_N = record_M(N,:)'*dt;
else
    for N= 1:i
        CI_temp = (record_M(N,:)'*dt); %CI is 3X1 matrix  
        CI = CI + CI_temp;
    end
   Z_N = [Z_N;CI]; %Z_N is 3Nx1 matrix, in here,CI should be transport.
end

GF = [0,0,0];
CF = zeros(3,6);
CI = [0,0,0]';

% using least square to estimate X
 X_hat = inv(Psi_N'*Psi_N)*Psi_N'*Z_N;

 J_hat = X_hat(1:6,:);
 
 r_CM_hat = X_hat(7:9,:);

 record_MOI(i,:) = X_hat(1:3,:);
 record_POI(i,:) = X_hat(4:6,:);
 record_CM(i,:) = r_CM_hat;


 end

%% Plot


% Omega plot

figure;

tiledlayout(3, 1);

ax1 = nexttile;
plot(ax1, ...  
          t, record_Omega(1:length(record_Omega),1),'-'...
          );
%title("Roll", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("\omega_{sys,x}", FontSize=13);
legend("x");

ax2 = nexttile;
plot(ax2, ...  
          t, record_Omega(1:length(record_Omega),2),'-'...
          );
%title("Pitch", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("\omega_{sys,y}", FontSize=13);
legend("y");


ax3 = nexttile;

plot(ax3, ...  
          t, record_Omega(1:length(record_Omega),3),'-'...
          );
%title("Yaw", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("\omega_{sys,z}", FontSize=13);
legend("z");

sgtitle('Platform Angular Rates (rad/sec)', 'FontSize', 16);


% Attitude plot

figure;

tiledlayout(3, 1);

ax4 = nexttile;
plot(ax4, ...  
          t, record_R(1:length(record_R),1),'-'...
          );
%title("Roll", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("\phi", FontSize=13);
legend("x");

ax5 = nexttile;
plot(ax5, ...  
          t, record_R(1:length(record_R),2),'-'...
          );
%title("Pitch", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("\theta", FontSize=13);
legend("y");


ax6 = nexttile;

plot(ax6, ...  
          t, record_R(1:length(record_R),3),'-'...
          );
%title("Yaw", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("\psi", FontSize=13);
legend("z");

sgtitle('Platform Attitude (rad)', 'FontSize', 16);

% Mass times C.M offset(kg*m)

figure;

tiledlayout(3, 1);

ax7 = nexttile;
plot(ax7, ...  
          t, record_CM(1:length(record_CM),1),'-'...
          );
%title("Roll", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("$mr_{x}$", FontSize=13);
legend("x");

ax8 = nexttile;
plot(ax8, ...  
          t, record_CM(1:length(record_CM),2),'-'...
          );
%title("Pitch", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("$mr_{y}$", FontSize=13);
legend("y");


ax9 = nexttile;

plot(ax9, ...  
          t, record_CM(1:length(record_CM),3),'-'...
          );
%title("Yaw", FontSize=14);
xlabel("Time(10ms)", FontSize=13);
ylabel("$mr_{z}$", FontSize=13);
legend("z");

sgtitle('Mass times C.M offset(kg*m)', 'FontSize', 16);

% M.O.I
figure;

tiledlayout(3, 1);

ax10 = nexttile;
plot(ax10, ...  
          t, record_MOI(1:length(record_MOI),1),'-'...
          );

xlabel("Time(10ms)", FontSize=13);
ylabel("$I_{x}$", FontSize=13);
legend("x");

ax11 = nexttile;
plot(ax11, ...  
          t, record_MOI(1:length(record_MOI),2),'-'...
          );

xlabel("Time(10ms)", FontSize=13);
ylabel("$I_{y}$", FontSize=13);
legend("y");


ax12 = nexttile;

plot(ax12, ...  
          t, record_MOI(1:length(record_MOI),3),'-'...
          );

xlabel("Time(10ms)", FontSize=13);
ylabel("$I_{z}$", FontSize=13);
legend("z");

sgtitle('Moment of Inertia(kg*m^2)', 'FontSize', 16);






% P.O.I
figure;

tiledlayout(3, 1);

ax13 = nexttile;
plot(ax13, ...  
          t, record_POI(1:length(record_POI),1),'-'...
          );

xlabel("Time(10ms)", FontSize=13);
ylabel("$I_{xy}$", FontSize=13);
legend("x");

ax14 = nexttile;
plot(ax14, ...  
          t, record_POI(1:length(record_POI),2),'-'...
          );

xlabel("Time(10ms)", FontSize=13);
ylabel("$I_{xz}$", FontSize=13);
legend("y");


ax15 = nexttile;

plot(ax15, ...  
          t, record_POI(1:length(record_POI),3),'-'...
          );

xlabel("Time(10ms)", FontSize=13);
ylabel("$I_{yz}$", FontSize=13);
legend("z");

sgtitle('Products of Inertia(kg*m^2)', 'FontSize', 16);













