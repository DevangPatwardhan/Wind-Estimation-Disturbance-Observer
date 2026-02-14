function v_rel= draginversion(u)


d_hat = [u(1); u(2); u(3)]; 
global rho Cd A 

eps_d = 1e-6;                 % small threshold to avoid sqrt(0)
d_mag = max(abs(d_hat), eps_d);

v_rel = sign(d_hat) .* sqrt( 2*d_mag ./ (rho*Cd*A) );

