function [ew] = Omega_Error(O,Od,R,Rd)
%OMEGA_ERROR Summary of this function goes here
%   Detailed explanation goes here
ew = [0;0;0];
ew = O - (R'*Rd*Od');
%ew = O -Od';
end

