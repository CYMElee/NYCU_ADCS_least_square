function [Od,Odd] = Get_Omega_Desire(rd,Od_raw,Odd_raw,i,dt)
%R is current Attitude
% Od_raw is desire angular velocity but not compare with Euler angle
% principle
Od = [0;0;0];
Odd = [0;0;0];

convert_martix = [1,0,-sin(rd(2));...
                  0,cos(rd(1)),sin(rd(1))*cos(rd(2));...
                  0,-sin(rd(1)),cos(rd(1))*cos(rd(2))];
Od = convert_martix*Od_raw';
Odd =convert_martix*Odd_raw';
end

