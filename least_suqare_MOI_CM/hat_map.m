function [Skew_symmetric_matrix] = hat_map(M)
Skew_symmetric_matrix = zeros(3);

Skew_symmetric_matrix = [0,-M(3),M(2);...
                         M(3),0,-M(1);...
                         -M(2),M(1),0];
end

