function [Skew_symmetric_matrix] = hat_map(T)
Skew_symmetric_matrix = zeros(3);

Skew_symmetric_matrix = [0,-T(3),T(2);...
                         T(3),0,-T(1);...
                         -T(2),T(1),0];
end

