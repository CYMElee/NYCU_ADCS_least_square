function [outputArg1] = Euler_matrix(inputArg1)
outputArg1 = eye(3);
outputArg1 = [1 ,sin(inputArg1(1))*tan(inputArg1(2)),cos(inputArg1(1))*tan(inputArg1(2)); ...
              0,cos(inputArg1(1)),-sin(inputArg1(1)) ; ...
              0,sin(inputArg1(1))*sec(inputArg1(2)) ,cos(inputArg1(1))*sec(inputArg1(2)) ];

end

