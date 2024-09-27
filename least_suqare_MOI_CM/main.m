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

record_Omega = zeros(length(t),3);

record_MOI_TRUE = zeros(length(t),3);
record_MOI = zeros(length(t),3); % a length(t)*3 vector using store estimated MOI

record_POI_TRUE = zeros(length(t),3);
record_POI = zeros(length(t),3); % a length(t)*3 vector using store estimated products of inertia

record_CM_TRUE = zeros(length(t),3);
record_CM = zeros(length(t),3); % a length(t)*3 vector using store estimated C.M offset

% torque comp
record_TOR_TRUE = zeros(length(t),3);
record_TOR_CMD = zeros(length(t),3);


% Psi record
record_centripetal_force = zeros(3,6,length(t));

record_gravitational_acceleration = zeros(3,3,length(t));

% M record replace by 2 turn 
record_Coriolis = zeros(length(t),3);
record_m = zeros(length(t),3);


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
    
    % if Euler angle's x,y belong(-5,5) we assign torque directly
    a1 = sin(r(3)/2)*cos(r(2)/2)*cos(r(1)/2)-cos(r(3)/2)*sin(r(2)/2)*sin(r(1)/2);
    a2 = cos(r(3)/2)*sin(r(2)/2)*cos(r(1)/2)+sin(r(3)/2)*cos(r(2)/2)*sin(r(1)/2);
    a3 = cos(r(3)/2)*cos(r(2)/2)*sin(r(1)/2)-sin(r(3)/2)*sin(r(2)/2)*cos(r(1)/2);
    alpha = [a1; ...
             a2; ...   
             a3];
    alpha_4 = cos(r(3)/2)*cos(r(2)/2)*cos(r(1)/2)+sin(r(3)/2)*sin(r(2)/2)*sin(r(1)/2);
    [alpha_hat] = hat_map(alpha);

    %else
    if  abs(r(3)) < 0.085 && abs(r(2)) < 0.085
        M_p = 2*[sin((2*pi*i*dt)/5);cos((2*pi*i*dt)/5);-sin((2*pi*i*dt)/5)];
       
    else
        M_p = -0.5*(((alpha_hat)+alpha_4*eye(3))*Gp + Gamma*(1-alpha_4)*eye(3))*alpha-Gr*omega_ab;
    end
    %record_M(i,:) = M_p';
    record_TOR_CMD(i,:)= M_p';
    M = [M_p;0];
    %% using A.B desire M to get R.W generate M
    omega_dot_mo = -(inv(H_w)/J_RW_testbed)*M;

    for m=1:4
        if (J_RW_testbed*omega_dot_mo(m)>1.7794)
            omega_dot_mo(m)= 1.7794/J_RW_testbed;
        end
        if (J_RW_testbed*omega_dot_mo(m)<-1.7794)
            omega_dot_mo(m)= -1.7794/J_RW_testbed;
        end
    end


    omega_mo = omega_mo_prev + omega_dot_mo*dt;
    %constrain the r.w rpm
    for m=1:4
        if (omega_mo(m)>=592)
            omega_mo(m)= 592;
        end
        if (omega_mo(m)<-592)
            omega_mo(m)= -592;
        end
    end
    record_m(i,:) = -(A_w*J_RW_testbed*(omega_mo-omega_mo_prev ));
    omega_mo_prev = omega_mo; 

    %% system dynamics
    % using R.C as momentum exchange devices.
    record_M(i,:) = - A_w*J_RW_testbed*omega_dot_mo -cross(omega_ab_prev,A_w*J_RW_testbed*omega_mo);
    record_Coriolis(i,:) = -cross(omega_ab_prev,A_w*J_RW_testbed*omega_mo);


    record_TOR_TRUE(i,:) =- A_w*J_RW_testbed*omega_dot_mo -cross(omega_ab_prev,A_w*J_RW_testbed*omega_mo);
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
 Psi_N = [omega_N-omega_init+CF ,GF ];




 %/*************************/%
 %/            Z_N          /%
 %/*************************/%

if i == 1
    Z_N = (record_m(1,:)'+(record_Coriolis(1,:)'*dt));
else
    for N= 1:i
        CI_temp = (record_m(N,:)'+(record_Coriolis(N,:)'*dt)); %CI is 3X1 matrix  
        CI = CI + CI_temp;
    end
    Z_N =CI ; %Z_N is 3Nx1 matrix, in here,CI should be transport.
end





GF = [0,0,0];
CF = zeros(3,6);
CI = [0,0,0]';

% using least square to estimate X

if i < 10
 X_hat = inv(Psi_N'*Psi_N)*Psi_N'*Z_N;
 X_hat_prev = X_hat;
 x_init = X_hat;
 z_init = Z_N;
 psi_init = Psi_N;
else
  % Gain Calculation   
  K_N =  P_N_prev*Psi_N'*inv(eye(3)+Psi_N*P_N_prev*Psi_N');
  % Covariance Update
  X_hat =  X_hat_prev + K_N*(Z_N-Psi_N*X_hat_prev);
  X_hat_prev = X_hat;
  P_N = (P_N_prev - (K_N*Psi_N*P_N_prev));
  P_N_prev = P_N;
end
 J_hat = X_hat(1:6,:);
 
 r_CM_hat = X_hat(7:9,:);

 record_MOI(i,:) = X_hat(1:3,:);
 record_POI(i,:) = X_hat(4:6,:);
 record_CM(i,:) = r_CM_hat;
 record_MOI_TRUE(i,:) =J_MOI;
 record_POI_TRUE(i,:) = J_POI;
 record_CM_TRUE(i,:) = r_CMXM;

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
xlabel("Time(10ms)", 'FontSize',13);
ylabel("\omega_{sys,x}",'FontSize',13);
legend("x");

ax2 = nexttile;
plot(ax2, ...  
          t, record_Omega(1:length(record_Omega),2),'-'...
          );
%title("Pitch", FontSize=14);
xlabel("Time(10ms)",'FontSize',13);
ylabel("\omega_{sys,y}", 'FontSize',13);
legend("y");


ax3 = nexttile;

plot(ax3, ...  
          t, record_Omega(1:length(record_Omega),3),'-'...
          );
%title("Yaw", FontSize=14);
xlabel("Time(10ms)", 'FontSize',13);
ylabel("\omega_{sys,z}", 'FontSize',13);
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
xlabel("Time(10ms)", 'FontSize',13);
ylabel("\psi", 'FontSize',13);
legend("x");

ax5 = nexttile;
plot(ax5, ...  
          t, record_R(1:length(record_R),2),'-'...
          );
%title("Pitch", FontSize=14);
xlabel("Time(10ms)", 'FontSize',13);
ylabel("\theta", 'FontSize',13);
legend("y");


ax6 = nexttile;

plot(ax6, ...  
          t, record_R(1:length(record_R),3),'-'...
          );
%title("Yaw", FontSize=14);
xlabel("Time(10ms)", 'FontSize',13);
ylabel("\phi", 'FontSize',13);
legend("z");

sgtitle('Platform Attitude (rad)', 'FontSize', 16);

% Mass times C.M offset(kg*m)

figure;

tiledlayout(3, 1);

ax7 = nexttile;
plot(ax7, t, record_CM(1:length(record_CM),1),'--',t, record_CM_TRUE(1:length(record_CM_TRUE),1),'-');
%title("Roll", FontSize=14);
xlabel("Time(10ms)", 'FontSize',13);
ylabel("mr_{x}", 'FontSize',13);
legend({'$\hat{mr_{x}}$', '$mr_{x}$'}, 'Interpreter', 'latex');

ax8 = nexttile;
plot(ax8, t, record_CM(1:length(record_CM),2),'--',t, record_CM_TRUE(1:length(record_CM_TRUE),2),'-');
%title("Pitch", FontSize=14);
xlabel("Time(10ms)", 'FontSize',13);
ylabel("mr_{y}", 'FontSize',13);
legend({'$\hat{mr_{y}}$', '$mr_{y}$'}, 'Interpreter', 'latex');


ax9 = nexttile;

plot(ax9, t, record_CM(1:length(record_CM),3),'--',t, record_CM_TRUE(1:length(record_CM_TRUE),3),'-');
%title("Yaw", FontSize=14);
xlabel("Time(10ms)", 'FontSize',13);
ylabel("mr_{z}", 'FontSize',13);
legend({'$\hat{mr_{z}}$', '$mr_{z}$'}, 'Interpreter', 'latex');

sgtitle('Mass times C.M offset(kg*m)', 'FontSize', 16);

% M.O.I
figure;

tiledlayout(3, 1);

ax10 = nexttile;
plot(ax10, t, record_MOI(1:length(record_MOI),1),'--',t, record_MOI_TRUE(1:length(record_MOI_TRUE),1),'-');

xlabel("Time(10ms)", 'FontSize',13);
ylabel("J_{x}", 'FontSize',13);
legend({'$\hat{J_{x}}$', '$J_{x}$'}, 'Interpreter', 'latex');

ax11 = nexttile;
plot(ax11, t, record_MOI(1:length(record_MOI),2),'--',t, record_MOI_TRUE(1:length(record_MOI_TRUE),2),'-');
xlabel("Time(10ms)", 'FontSize',13);
ylabel("J_{y}", 'FontSize',13);
legend({'$\hat{J_{y}}$', '$J_{y}$'}, 'Interpreter', 'latex');


ax12 = nexttile;

plot(ax12, t, record_MOI(1:length(record_MOI),3),'--',t, record_MOI_TRUE(1:length(record_MOI_TRUE),3),'-');

xlabel("Time(10ms)", 'FontSize',13);
ylabel("J_{z}", 'FontSize',13);
legend({'$\hat{J_{z}}$', '$J_{z}$'}, 'Interpreter', 'latex');

sgtitle('Moment of Inertia(kg*m^2)', 'FontSize', 16);






% P.O.I
figure;

tiledlayout(3, 1);

ax13 = nexttile;
plot(ax13, t, record_POI(1:length(record_POI),1),'--',t, record_POI_TRUE(1:length(record_POI_TRUE),1),'-');
xlabel("Time(10ms)", 'FontSize',13);
ylabel("J_{xy}", 'FontSize',13);
legend({'$\hat{J_{xy}}$', '$J_{xy}$'}, 'Interpreter', 'latex');

ax14 = nexttile;
plot(ax14, t, record_POI(1:length(record_POI),2),'--',t, record_POI_TRUE(1:length(record_POI_TRUE),2),'-');

xlabel("Time(10ms)",'FontSize',13);
ylabel("J_{xz}", 'FontSize',13);
legend({'$\hat{J_{xz}}$', '$J_{xz}$'}, 'Interpreter', 'latex');


ax15 = nexttile;

plot(ax15, t, record_POI(1:length(record_POI),3),'--',t, record_POI_TRUE(1:length(record_POI_TRUE),3),'-');

xlabel("Time(10ms)", 'FontSize',13);
ylabel("J_{yz}", 'FontSize',13);
legend({'$\hat{J_{yz}}$', '$J_{yz}$'}, 'Interpreter', 'latex');

sgtitle('Products of Inertia(kg*m^2)', 'FontSize', 16);



% Torque vs Torque_cmd
figure;
ax16 = nexttile;

plot(ax16, t, record_TOR_TRUE(1:length(record_TOR_TRUE),1),'--',t, record_TOR_CMD(1:length(record_TOR_CMD),1),'-');

xlabel("Time(10ms)", 'FontSize',13);
ylabel("Mx", 'FontSize',13);
legend({'$Mx\_truth$', '$Mx\_command$'}, 'Interpreter', 'latex');
sgtitle('Torque Generate by R.W(x)', 'FontSize', 16);

figure;

ax17 = nexttile;
plot(ax17, t, record_TOR_TRUE(1:length(record_TOR_TRUE),2),'--',t, record_TOR_CMD(1:length(record_TOR_CMD),2),'-');
xlabel("Time(10ms)", 'FontSize',13);
ylabel("My", 'FontSize',13);
legend({'$My\_truth$', '$My\_command$'}, 'Interpreter', 'latex');
sgtitle('Torque Generate by R.W(y)', 'FontSize', 16);

figure;
ax18 = nexttile;
plot(ax18, t, record_TOR_TRUE(1:length(record_TOR_TRUE),3),'--',t, record_TOR_CMD(1:length(record_TOR_CMD),3),'-');
xlabel("Time(10ms)", 'FontSize',13);
ylabel("Mz", 'FontSize',13);
legend({'$Mz\_truth$', '$Mz\_command$'}, 'Interpreter', 'latex');
sgtitle('Torque Generate by R.W(z)', 'FontSize', 16);


%     %error


Percent_error_MOI_x =100*((record_MOI_TRUE(length(record_MOI_TRUE),1)-record_MOI(length(record_MOI),1))/record_MOI_TRUE(length(record_MOI_TRUE),1));
Percent_error_MOI_y =100*((record_MOI_TRUE(length(record_MOI_TRUE),2)-record_MOI(length(record_MOI),2))/record_MOI_TRUE(length(record_MOI_TRUE),2));
Percent_error_MOI_z =100*((record_MOI_TRUE(length(record_MOI_TRUE),3)-record_MOI(length(record_MOI),3))/record_MOI_TRUE(length(record_MOI_TRUE),3));


Percent_error_POI_x = 100*((record_POI_TRUE(length(record_POI_TRUE),1)-record_POI(length(record_POI),1))/record_POI_TRUE(length(record_POI_TRUE),1));
Percent_error_POI_y = 100*((record_POI_TRUE(length(record_POI_TRUE),2)-record_POI(length(record_POI),2))/record_POI_TRUE(length(record_POI_TRUE),2));
Percent_error_POI_z = 100*((record_POI_TRUE(length(record_POI_TRUE),3)-record_POI(length(record_POI),3))/record_POI_TRUE(length(record_POI_TRUE),3));


Percent_error_CM_x = 100*((record_CM_TRUE(length(record_CM_TRUE),1)-record_CM(length(record_CM),1))/record_CM_TRUE(length(record_CM),1));
Percent_error_CM_y = 100*((record_CM_TRUE(length(record_CM_TRUE),2)-record_CM(length(record_CM),2))/record_CM_TRUE(length(record_CM),2));
Percent_error_CM_z = 100*((record_CM_TRUE(length(record_CM_TRUE),3)-record_CM(length(record_CM),3))/record_CM_TRUE(length(record_CM),3));









