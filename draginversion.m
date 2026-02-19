function v_rel= draginversion(u)


d_hat = [u(1); u(2); u(3)]; 
% Force dimension
d_hat = reshape(d_hat,3,1);

v_rel = zeros(3,1);   % preallocate fixed size

rho = 1.225;
Cd  = 1.0;
A   = 0.1;

k = 0.5*rho*Cd*A;

d_norm = sqrt(d_hat.'*d_hat);

if d_norm > 1e-6
    v_mag = sqrt(d_norm/k);
    v_rel(:) = v_mag * (d_hat/d_norm);
end
end