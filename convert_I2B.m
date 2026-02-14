function vector_body = convert_I2B(vector_inertial, state)
% CONVERT_I2B Converts a vector from inertial frame to body frame using Euler angles.
%
% INPUTS:
%   vector_inertial - (3x1 vector) Vector in the inertial frame (e.g., wind velocity)
%   state           - (nx1 vector) State vector containing orientation angles:
%                     state(4) = phi (roll), state(5) = theta (pitch), state(6) = psi (yaw)
%
% OUTPUT:
%   vector_body     - (3x1 vector) Vector transformed to body frame

    % Extract Euler angles from the state vector
    phi = state(4);    % Roll angle
    theta = state(5);  % Pitch angle
    psi = state(6);    % Yaw angle

    % Precompute sines and cosines of the Euler angles
    sp = sin(phi);     % sin(phi)
    cp = cos(phi);     % cos(phi)
    st = sin(theta);   % sin(theta)
    ct = cos(theta);   % cos(theta)
    ss = sin(psi);     % sin(psi)
    cs = cos(psi);     % cos(psi)

    % Rotation matrix from body frame to inertial frame (Z-Y-X Euler sequence)
    R_B2I = [ ct*cs,  sp*st*cs - cp*ss,  cp*st*cs + sp*ss;
              ct*ss,  sp*st*ss - cp*cs,  cp*st*ss - sp*cs;
             -st,     sp*ct,            cp*ct ];

    % Transpose the rotation matrix to get inertial-to-body transformation
    R_I2B = R_B2I';

    % Transform the vector from inertial frame to body frame
    vector_body = R_I2B * vector_inertial;
end
