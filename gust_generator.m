function W_I = gust_generator(t)

% Inputs:
%   t     - time (seconds)
%   mode  - gust type:
%           1 = step
%           2 = ramp
%           3 = 1-cosine
%           4 = ramp (fixed direction x-z plane)
%           5 = overlay (ramp + cosine)
%
% Output:
%   W_I   - 3x1 wind vector in INERTIAL frame

mode = 5;

W_I = [0;0;0];

switch mode

    %% 1. Step gust (axial z-direction)
    case 1
        t0 = 2;
        W0 = [0;0;5];
        if t >= t0
            W_I = W0;
        end

    %% 2. Ramp gust (axial z-direction)
    case 2
        t0 = 2; 
        Tr = 3;
        Wmax = [0;0;5];

        if t >= t0 && t <= t0+Tr
            W_I = Wmax*(t-t0)/Tr;
        elseif t > t0+Tr
            W_I = Wmax;
        end

    %% 3. 1-Cosine gust (axial z-direction)
    case 3
        t0 = 3; 
        Tg = 30;
        Wmax = [0;0;3];

        if t >= t0 && t <= t0+Tg
            W_I = 0.5*Wmax*(1 - cos(pi*(t-t0)/Tg));
        elseif t > t0+Tg
            W_I = Wmax;
        end

    %% 4. Ramp with fixed direction in x-z plane
    case 4
        t0 = 2; 
        Tr = 3;
        theta = pi/6;      % fixed direction
        dir = [sin(theta); 0; cos(theta)];
        Wmax = 4*dir;

        if t >= t0 && t <= t0+Tr
            W_I = Wmax*(t-t0)/Tr;
        elseif t > t0+Tr
            W_I = Wmax;
        end

    %% 5. Overlay: ramp + cosine (different directions)
    case 5
        % Directions
        theta1 = pi/4;
        theta2 = -pi/6;

        dir1 = [sin(theta1); 0; cos(theta1)];
        dir2 = [sin(theta2); 0; cos(theta2)];

        % --- Ramp component ---
        t0 = 2; 
        Tr = 3;
        W1 = [0;0;0];

        if t >= t0 && t <= t0+Tr
            W1 = 3*dir1*(t-t0)/Tr;
        elseif t > t0+Tr
            W1 = 3*dir1;
        end

        % --- Cosine component ---
        t1 = 6; 
        Tg = 2;
        W2 = [0;0;0];

        if t >= t1 && t <= t1+Tg
            W2 = 2*dir2*0.5*(1 - cos(pi*(t-t1)/Tg));
        elseif t > t1+Tg
            W2 = 2*dir2;
        end

        W_I = W1 + W2;


end

end
