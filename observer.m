function d_hat_dot = observer(u)

d_hat_dot = zeros(3,1);

% Gains (scalar)
L1 = 3;
L2 = 3;
L3 = 5;

% Error signals (scalars)
e1 = u(1);
e2 = u(2);
e3 = u(3);

% Observer equations 
d_hat_dot(1) = L1 * e1;
d_hat_dot(2) = L2 * e2;
d_hat_dot(3) = L3 * e3;

end