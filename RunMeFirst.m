clear; clc; close all

global ...
    m g ...
    Ixx Iyy Izz Ixy Ixz Iyz ...
    inertial_matrix inv_I ...
    n_rotor rotor_l n_states ...
    rho Cd A ...
    use_bemt ...
    L_c T_c M_c F_c ...
    m_nd g_nd ...
    Ixx_nd Iyy_nd Izz_nd Ixy_nd Ixz_nd Iyz_nd ...
    inv_I_nd rotor_l_nd rho_nd A_nd ...
    R A_rotor ...
    R_nd A_rotor_nd Ct Cq ...
    n_blades n_segments theta_r ...
    cl_alpha cd_0 aoa_0 sigma r_cout ...
    c_nd c r_span_nd ...
    fwd_vel


m = 12;           % kg
g = 9.8;

Ixx = 0.6430;
Iyy = 0.6820;
Izz = 0.8320;
Ixy = 0;
Ixz = 0;
Iyz = 0;

inertial_matrix = [ ...
     Ixx  -Ixy  -Ixz; ...
    -Ixy   Iyy  -Iyz; ...
    -Ixz  -Iyz   Izz ];

inv_I = inv(inertial_matrix);

n_rotor  = 6;
rotor_l  = 0.6;
n_states = 12;

rho = 1.2256;
Cd  = 2;
A   = [0.3; 0.3; 0.3];

use_bemt = 0;

L_c = rotor_l;
T_c = sqrt(L_c/g);
M_c = m;
F_c = M_c*L_c/(T_c^2);

% Nondimensional aircraft params
m_nd   = m/M_c;
g_nd   = g*(T_c^2/L_c);

Ixx_nd = Ixx/(M_c*L_c^2);
Iyy_nd = Iyy/(M_c*L_c^2);
Izz_nd = Izz/(M_c*L_c^2);
Ixy_nd = Ixy/(M_c*L_c^2);
Ixz_nd = Ixz/(M_c*L_c^2);
Iyz_nd = Iyz/(M_c*L_c^2);

inv_I_nd = inv_I*(M_c*L_c^2);

rotor_l_nd = rotor_l/L_c;
rho_nd     = rho*(L_c^3/M_c);
A_nd       = A/(L_c^2);


R_nd       = 0.2286/L_c;
A_rotor_nd = pi*R_nd^2;

R = 0.2286;                % rotor radius (m)
A_rotor = pi*R^2;

Ct = 0.008;
Cq = Ct*sqrt(Ct/2);

n_blades   = 2;
n_segments = 10;

theta_r  = 12*pi/180;
cl_alpha = 2*pi;
cd_0     = 0.01;
aoa_0    = 0;

sigma  = 0.042;
r_cout = 0.15;

c_nd = sigma*pi*R_nd/n_blades;
c = sigma*pi*R/n_blades;


r_span_nd = zeros(n_segments,1);
dy  = R_nd*(1-r_cout)/n_segments;
fac = r_cout*R_nd + 0.5*dy;

for k = 1:n_segments
    r_span_nd(k) = fac + (k-1)*dy;
end
fwd_vel = 5;


