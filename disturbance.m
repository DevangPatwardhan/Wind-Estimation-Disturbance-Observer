function d_raw = disturbance(u)

global m 

F = [u(1); u(2); u(3)];
v_dot = [u(4); u(5); u(6)];


d_raw = m*v_dot - F;

