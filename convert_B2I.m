function vector_inertial = convert_B2I(vector_body, state)
% CONVERT_B2I Converts a vector from body frame to inertial frame using Euler angles.
%
% INPUTS:
%   vector_body    - (3x1 vector) Vector in body frame (e.g., force, velocity)
%   state          - (nx1 vector) State vector containing orientation angles:
%                    state(4) = phi (roll), state(5) = theta (pitch), state(6) = psi (yaw)
%
% OUTPUT:
%   vector_inertial - (3x1 vector) Vector transformed to inertial frame

    % Extract Euler angles from the state vector
    phi = state(4);    % Roll angle
    theta = state(5);  % Pitch angle
    psi = state(6);    % Yaw angle

    % Precompute sines and cosines of the Euler angles
    sp = sin(phi);     
    cp = cos(phi);     
    st = sin(theta);   
    ct = cos(theta);   
    ss = sin(psi);     
    cs = cos(psi);     

    % Rotation matrix from body frame to inertial frame (Z-Y-X Euler sequence)
    R_B2I = [ ct*cs,  sp*st*cs - cp*ss,  cp*st*cs + sp*ss;
              ct*ss,  sp*st*ss - cp*cs,  cp*st*ss - sp*cs;
             -st,     sp*ct,            cp*ct ];

    % Transform the vector from body frame to inertial frame
    vector_inertial = R_B2I * vector_body;
end
