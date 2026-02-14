function FT = dynamics(u)
% DYNAMICS
% Unified force–moment model for a hexarotor
% Models supported:
%   1) Quadratic
%   2) BEMT (Uniform inflow)
%   3) Pitt–Peters
%
% INPUT:
%   csvfile : CSV file containing
%     t x y z phi theta psi u v w p q r u1 u2 u3 u4 u5 u6
%
% OUTPUT:
%   FT : [Fx Fy Fz Tx Ty Tz]'  (6x1) in BODY frame

% Load data


x     = u(1);
y     = u(2);
z     = u(3);
phi   = u(4);
theta = u(5);
psi   = u(6);
ub    = u(7);
vb    = u(8);
wb    = u(9);
p     = u(10);
q     = u(11);
r     = u(12);
input = u(13:18);   % rotor speeds


% State vectors
state = zeros(12,1);
state(1:3)   = [x; y; z];
state(4:6)   = [phi; theta; psi];
state(7:9)   = [ub; vb; wb];
state(10:12) = [p; q; q];   % <-- FIXED below
vel_b = state(7:9);


% Load global params

global n_rotor rotor_l 

%% ------------------------------------------------------------------------
% Rotor configuration
%% ------------------------------------------------------------------------
rotor_angles = (0:n_rotor-1) * (2*pi/n_rotor);
rotor_pos = [ ...
    sin(rotor_angles);
    cos(rotor_angles);
    zeros(1,n_rotor)
] * rotor_l;

rotor_spin = ones(n_rotor,1);
rotor_spin(2:2:end) = -1;



%% ------------------------------------------------------------------------
% Rotor force & torque
%% ------------------------------------------------------------------------
rotor_T = zeros(n_rotor,3);
rotor_Q = zeros(n_rotor,3);

global rho


global R_nd A_rotor_nd Ct Cq
R  = R_nd;
A  = A_rotor_nd;

% QUADRATIC MODEL

T = Ct * rho * A .* (input.^2) * R^2;
Q = Cq * rho * A .* (input.^2) * R^3;

rotor_T(:,3) = T;
rotor_Q(:,3) = Q;



% Total force & moment at CG
T_b = sum(rotor_T,1)';

M_b = zeros(3,1);
for i = 1:n_rotor
    M_b = M_b + cross(rotor_pos(:,i), rotor_T(i,:)') ...
              + rotor_spin(i)*rotor_Q(i,:)';
end
global m g
% Weight
W_I = [0;0;-m*g];
W_b = convert_I2B(W_I, state);

global Cd 
% Drag
D_b = -0.5 * rho *Cd * A .* abs(vel_b) .* vel_b;

% Final output
Fcg_b = T_b + W_b + D_b;

FT = [Fcg_b; M_b];

end
