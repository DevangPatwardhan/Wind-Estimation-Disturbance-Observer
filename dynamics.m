function FT = dynamics(u)

% Fully dimensional rotor force–moment model
% Returns [Fx Fy Fz Tx Ty Tz]' in BODY frame (SI units)


x = u(1);
y = u(2);
z = u(3);


phi   = u(4);
theta = u(5);
psi   = u(6);

ub = u(7);
vb = u(8);
wb = u(9);

p = u(10);
q = u(11);
r = u(12);

Omega = u(13:18);   % rad/s (DIMENSIONAL)

R_ib = eul2rotm([psi, theta, phi]);  
vel_b   = [ub; vb; wb];
omega_b = [p; q; r];


global n_rotor rotor_l
global rho R A_rotor
global Ct Cq
global use_bemt
global n_blades n_segments theta_r
global cl_alpha cd_0 aoa_0
global c r_span

% Rotor Geometry

angles = (0:n_rotor-1)*(2*pi/n_rotor);

rotor_pos = [ ...
    rotor_l*cos(angles);
    rotor_l*sin(angles);
    zeros(1,n_rotor)];

spin = ones(n_rotor,1);
spin(2:2:end) = -1;


% Initialize totals
F_total = zeros(3,1);
M_total = zeros(3,1);

% LOOP OVER ROTORS

for i = 1:n_rotor
    
    pos_i   = rotor_pos(:,i);
    Omega_i = Omega(i);
    
    % Velocity at rotor disk (BODY frame)
    V_disk = vel_b + cross(omega_b,pos_i);
    Vz = V_disk(3);
    
    % MODEL SELECTION

    if use_bemt == 0
        
        % ---- KW^2 MODEL ----
        
        T = Ct * rho * A_rotor * (Omega_i * R)^2;
        Q = Cq * rho * A_rotor * (Omega_i * R)^2 * R;
        
    elseif use_bemt == 1
        
        % ---- Uniform Inflow Momentum Theory ----
        
        vi = 0;
        
        for k = 1:25
            
            T_mom = 2*rho*A_rotor*vi*(Vz + vi);
            
            lambda = (Vz + vi)/(Omega_i*R + 1e-6);
            CT_be  = Ct*(1 - lambda);
            
            T_be = CT_be * rho * A_rotor * (Omega_i*R)^2;
            
            dT = 2*rho*A_rotor*(Vz + 2*vi);
            vi = vi - (T_mom - T_be)/(dT + 1e-6);
        end
        
        T = 2*rho*A_rotor*vi*(Vz + vi);
        Q = T * R * sqrt(abs(T)/(2*rho*A_rotor + 1e-6));
        
    else
        
        % ---- Radial BEMT ----
        
        T = 0;
        Q = 0;
        
        dr = r_span(2) - r_span(1);
        
        for k = 1:n_segments
            
            r_loc = r_span(k);
            
            Ut = Omega_i * r_loc;
            Ua = Vz;
            
            phi_i = atan2(Ua, Ut);
            
            alpha = theta_r - phi_i - aoa_0;
            
            Cl = cl_alpha * alpha;
            Cd = cd_0;
            
            Vrel2 = Ut^2 + Ua^2;
            
            dL = 0.5*rho*Vrel2*c*Cl*dr;
            dD = 0.5*rho*Vrel2*c*Cd*dr;
            
            dT = n_blades*( dL*cos(phi_i) - dD*sin(phi_i) );
            dQ = n_blades*( r_loc*(dL*sin(phi_i) + dD*cos(phi_i)) );
            
            T = T + dT;
            Q = Q + dQ;
        end
    end
    

    % Force & Moment Contribution

    
    F_i = [0; 0; -T];          % thrust along -Z body
    % M_offset = cross(pos_i, F_i);
    % M_react  = [0; 0; spin(i)*Q];
    
    F_total = F_total + F_i;
    M_total = M_total; %+ M_offset + M_react;
end

% Output


FT = [F_total; M_total];

end
