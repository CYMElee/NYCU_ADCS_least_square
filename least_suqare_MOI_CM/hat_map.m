function [Skew_symmetric_matrix] = hat_map(M)
Skew-symmetric_matrix = zeros(3);

Skew-symmetric_matrix = [0,-M(3,1),M(2,1);...
                         M(3),0,-M(1,1);...
                         -M(2),M(1),0];

end

