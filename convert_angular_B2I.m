function vector_inertial = convert_angular_B2I(vector_body, state)
% CONVERT_ANGULAR_B2I Converts angular velocity from body frame to inertial frame rates.
%
% INPUTS:
%   vector_body    - (3x1 vector) Angular velocity in body frame [p; q; r]
%   state          - (nx1 vector) State vector containing orientation
%                    angles: state(4) = phi (roll), state(5) = theta (pitch)
%
% OUTPUT:
%   vector_inertial - (3x1 vector) [phi_dot; theta_dot; psi_dot] in inertial frame

    % Extract roll (phi) and pitch (theta) angles from the state vector
    phi = state(4); 
    theta = state(5); 

    % Compute sines and cosines for transformation matrix
    sp = sin(phi);   
    cp = cos(phi);   
    tt = sin(theta); 
    ct = cos(theta); 

    % Transformation matrix from body angular rates [p; q; r] 
    % to inertial angular rates [phi_dot; theta_dot; psi_dot]
    T = [1, sp*tt,  cp*tt;
         0, cp,     -sp;
         0, sp/ct,  cp/ct];

    % Apply transformation
    vector_inertial = T * vector_body;
end
