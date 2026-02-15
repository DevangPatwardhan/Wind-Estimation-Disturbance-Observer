function Xdot = simulate_hexarotor(u)

global ...
    m g ...
    Ixx Iyy Izz Ixy Ixz Iyz ...
    inv_I ...
    n_rotor rotor_l ...
    rho Cd A ...
    R_nd A_rotor_nd Ct Cq
X      = u(1:18);
wind_I = u(19:21);


pos     = X(1:3);
phi     = X(4); 
theta   = X(5); 
psi     = X(6);

vel_b   = X(7:9);      % body velocity
omega   = X(10:12);    % p q r
omega_r = X(13:18);   % rotor speeds

p = omega(1); 
q = omega(2); 
r = omega(3);



R_ib = eul2rotm([psi theta phi]);   % body -> inertial
R_bi = R_ib';

vel_I = R_ib * vel_b;

wind_b = R_bi * wind_I;

V_rel_b = vel_b - wind_b;



flight_mode = 'climb';   

switch flight_mode
    case 'hover'
        pos_des = [0;0;10];
        vel_des = [0;0;0];
    case 'climb'
        pos_des = [0;0;50];
        vel_des = [0;0;5];
    case 'forward'
        pos_des = [100;0;0];
        vel_des = [5;0;0];
end





psi_des = 0;

% --- Position loop (outer loop)
Kp_pos = 1.5;
Kd_pos = 2.5;    % damping term

% --- Attitude loop (inner loop)
Kp_att_rollpitch = 8;
Kd_att_rollpitch = 4.5;

Ki_pos = 0.1;  % tune



pos_err = pos_des - pos;
vel_err = vel_des - vel_I;

% IMP commmanding control eq -- (do not change this, change the gains instead- Devang)
v_cmd = Kp_pos * pos_err;
v_max = 5;  % or whatever you want
v_cmd = max(min(v_cmd, v_max), -v_max);
a_cmd = Kd_pos * (v_cmd - vel_I);


F_des_I = m*(a_cmd + [0;0;g]);
T_total = dot(F_des_I, R_ib(:,3));

zb_des = F_des_I / (norm(F_des_I) + 1e-6);

% Desired heading direction
xc_des = [cos(psi_des); sin(psi_des); 0];

yb_des = cross(zb_des, xc_des);
yb_des = yb_des / (norm(yb_des) + 1e-6);

xb_des = cross(yb_des, zb_des);

R_des = [xb_des yb_des zb_des];

% Current rotation
R = R_ib;

% Rotation error (Lee et al. style)
e_R_mat = 0.5*(R_des' * R - R' * R_des);
e_R = [e_R_mat(3,2);
       e_R_mat(1,3);
       e_R_mat(2,1)];

e_omega = omega;   % desired body rates = 0

% Control law
Kp_att = diag([Kp_att_rollpitch Kp_att_rollpitch 2.0]);
Kd_att = diag([Kd_att_rollpitch Kd_att_rollpitch 1.5]);

tau_cmd = -Kp_att*e_R - Kd_att*e_omega;





L = rotor_l;

angles = deg2rad([0 60 120 180 240 300]);

x_i = L * cos(angles);
y_i = L * sin(angles);

M = zeros(4,6);

for i = 1:6
    M(1,i) = 1;              % thrust
    M(2,i) = y_i(i);         % tau_x
    M(3,i) = -x_i(i);        % tau_y
    
    if mod(i,2)==0
        M(4,i) = -Cq;
    else
        M(4,i) =  Cq;
    end
end








desired = [T_total; tau_cmd];

T_i = pinv(M) * desired;
T_i = max(T_i,0);



F_thrust_total = sum(T_i);

M_reaction_total = 0;

for i = 1:n_rotor
    
    Ti = T_i(i);
    
    % compute rotor speed from thrust (static model)
    omega_i = sqrt( Ti / (rho*A_rotor_nd*Ct*R_nd^2 + 1e-6) );
    
    Qi = rho*A_rotor_nd*(omega_i*R_nd)^2 * Cq;
    
    if mod(i,2)==0
        M_reaction_total = M_reaction_total - Qi;
    else
        M_reaction_total = M_reaction_total + Qi;
    end
end


F_thrust_b = [0;0;F_thrust_total];



F_drag = -0.5*rho*Cd*A .* abs(V_rel_b).*V_rel_b;


F_g_b = R_bi * [0;0;-m*g];



F_b = F_thrust_b + F_drag + F_g_b;



vel_dot = (1/m)*F_b - cross(omega, vel_b);



I = [ Ixx -Ixy -Ixz;
     -Ixy  Iyy -Iyz;
     -Ixz -Iyz  Izz ];

M_total = [tau_cmd(1);
           tau_cmd(2);
           tau_cmd(3) + M_reaction_total];

omega_dot = inv_I * ( M_total - cross(omega, I*omega) );



T_eul = [
    1  sin(phi)*tan(theta)  cos(phi)*tan(theta);
    0  cos(phi)            -sin(phi);
    0  sin(phi)/cos(theta)  cos(phi)/cos(theta)
];

eul_dot = T_eul * omega;



pos_dot = vel_I;



tau_motor = 0.03;
omega_cmd = sqrt(T_i./(rho*A_rotor_nd*Ct*R_nd^2 + 1e-6));
omega_r_dot = (omega_cmd - omega_r)/tau_motor;

int_pos_err = pos_err;

Xdot = [
    pos_dot;
    eul_dot;
    vel_dot;
    omega_dot;
    omega_r_dot;
];

end
