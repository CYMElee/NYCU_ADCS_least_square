function [outputArg1] = Angular_rate_matrix(inputArg1)
%ANGULAR_RATE_MATRIX Summary of this function goes here
outputArg1 = zeros(3,6);
outputArg1 = [inputArg1(1),0,0,inputArg1(2),inputArg1(3),0; ...
              0,inputArg1(2),0,inputArg1(1),0,inputArg1(3); ...
              0,0,inputArg1(3),0,inputArg1(1),inputArg1(3)];

end

