function [o_b] = omega_bar(omega_ab_prev,R,Rd,Od,Odd)
o_b =[0;0;0];
hat_map = [0,-omega_ab_prev(3),omega_ab_prev(2);...
           omega_ab_prev(3),0,-omega_ab_prev(1);...
           -omega_ab_prev(2),omega_ab_prev(1),0];
o_b_temp = (hat_map*R.'*Rd*Od')-(R.'*Rd*Odd.');

o_b = o_b_temp;
end

