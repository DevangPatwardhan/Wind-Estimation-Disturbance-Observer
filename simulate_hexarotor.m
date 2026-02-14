function [t, X] = simulate_hexarotor(wind_fun)

% wind_fun(t) -> returns wind vector in INERTIAL frame

global ...
    m g ...
    Ixx Iyy Izz Ixy Ixz Iyz ...
    inertial_matrix inv_I ...
    n_rotor rotor_l n_states ...
    rho Cd A ...
    R_nd A_rotor_nd Ct Cq ...
    inv_I_nd

%% -----------------------------
% SIMULATION SETTINGS
%% -----------------------------

dt = 0.002;
t0 = 0;
tf = 20;
t = t0:dt:tf;
N = length(t);

%% -----------------------------
% Flight Mode (DEFINE HERE)
%% -----------------------------

flight_mode = 'forward';   % 'hover','climb','forward'
inflow_mode = 'axial';     % 'axial','edgewise','none'

switch flight_mode
    case 'hover'
        V_des = [0;0;0];
    case 'climb'
        V_des = [0;0;3];
    case 'forward'
        V_des = [5;0;0];
end

%% -----------------------------
% State Definition
%% -----------------------------
% [x y z phi theta psi u v w p q r u1..u6]

X = zeros(12+n_rotor, N);

%% Initial rotor trim for hover
T_hover = m*g;
omega_trim = sqrt(T_hover/(n_rotor*Ct*rho*A_rotor_nd*R_nd^2));
X(13:end,1) = omega_trim;

%% Integral error
int_err = zeros(3,1);

%% Controller gains
Kp_pos = 3;
Ki_pos = 0.5;
Kp_att = 8;
Kd_att = 2;

%% -----------------------------
% MAIN LOOP
%% -----------------------------

for k = 1:N-1

    state = X(:,k);

    pos = state(1:3);
    phi = state(4); theta = state(5); psi = state(6);
    u = state(7); v = state(8); w = state(9);
    p = state(10); q = state(11); r = state(12);
    omega_r = state(13:end);

    vel_b = [u;v;w];
    omega = [p;q;r];

    %% Rotation matrix
    R = eul2rotm([psi theta phi]); % ZYX

    %% Wind
    wind_I = wind_fun(t(k));
    wind_b = R' * wind_I;

    vel_rel = vel_b - wind_b;

    %% Drag
    F_drag = -0.5*rho*Cd*A .* abs(vel_rel).*vel_rel;

    %% Velocity inertial
    vel_I = R * vel_b;

    %% Closed-loop velocity tracking
    vel_err = V_des - vel_I;
    int_err = int_err + vel_err*dt;

    a_cmd = Kp_pos*vel_err + Ki_pos*int_err;

    %% Desired thrust vector
    F_des_I = m*(a_cmd + [0;0;g]);
    T_mag = norm(F_des_I);

    %% Desired pitch from thrust direction
    theta_des = atan2(F_des_I(1),F_des_I(3));
    phi_des = -atan2(F_des_I(2),F_des_I(3));

    %% Attitude control torque
    tau = [
        -Kp_att*(phi-phi_des) - Kd_att*p;
        -Kp_att*(theta-theta_des) - Kd_att*q;
        -Kp_att*(psi) - Kd_att*r
    ];

    %% Rotor thrust with inflow model
    T_i = zeros(n_rotor,1);

    for i=1:n_rotor

        base_thrust = Ct*rho*A_rotor_nd*(omega_r(i)^2)*R_nd^2;

        switch inflow_mode

            case 'axial'
                ka = 0.5;
                T_i(i) = base_thrust - ka*vel_rel(3);

            case 'edgewise'
                ke = 0.3;
                V_edge = sqrt(vel_rel(1)^2 + vel_rel(2)^2);
                T_i(i) = base_thrust - ke*V_edge;

            otherwise
                T_i(i) = base_thrust;
        end
    end

    T_total = sum(T_i);
    F_thrust = [0;0;T_total];

    %% Gravity in body
    Fg_b = R' * [0;0;-m*g];

    %% Total forces
    F_b = F_thrust + F_drag + Fg_b;

    fx = F_b(1);
    fy = F_b(2);
    fz = F_b(3);

    %% -----------------------------
    % FULL BODY TRANSLATIONAL DYNAMICS
    %% -----------------------------

    u_dot = -(1/m)*(fx + r*v - q*w);
    v_dot = -(1/m)*(fy + p*w - r*u);
    w_dot = -(1/m)*(fz + q*u - p*v);

    %% -----------------------------
    % FULL BODY ROTATIONAL DYNAMICS
    %% -----------------------------

    mx = tau(1);
    my = tau(2);
    mz = tau(3);

    moment_vector = [
        mx + q*(Ixz*p + Iyz*q - Izz*r) - r*(Ixy*p - Iyy*q + Iyz*r);
        my - p*(Ixz*p + Iyz*q - Izz*r) + r*(Ixy*q - Ixx*p + Ixz*r);
        mz + p*(Ixy*p - Iyy*q + Iyz*r) - q*(Ixy*q - Ixx*p + Ixz*r)
    ];

    omega_dot = inv_I * moment_vector;

    %% Euler angle rates
    T_eul = [
        1  sin(phi)*tan(theta)  cos(phi)*tan(theta);
        0  cos(phi)            -sin(phi);
        0  sin(phi)/cos(theta)  cos(phi)/cos(theta)
    ];

    eul_dot = T_eul * omega;

    %% Position derivative
    pos_dot = vel_I;

    %% Motor first-order dynamics
    omega_cmd = sqrt(T_mag/(n_rotor*Ct*rho*A_rotor_nd*R_nd^2));
    omega_cmd = max(omega_cmd,0);
    tau_motor = 0.03;

    omega_r_dot = (omega_cmd - omega_r)/tau_motor;

    %% Assemble derivative
    Xdot = [
        pos_dot;
        eul_dot;
        u_dot;
        v_dot;
        w_dot;
        omega_dot;
        omega_r_dot
    ];

    %% RK4 Integration
    X(:,k+1) = X(:,k) + dt*Xdot;

end

end
