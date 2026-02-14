function W_I = gust_generator(t)

% Returns wind vector in INERTIAL frame
% Modify "mode" to switch gust types

mode = 5;   % 1=step, 2=ramp, 3=1cos, 4=random ramp, 5=overlay, 6=noisy overlay

W_I = [0;0;0];

switch mode

    %% 1. Step gust (axial, z-direction)
    case 1
        t0 = 2;
        W0 = [0;0;5];
        if t >= t0
            W_I = W0;
        end

    %% 2. Ramp gust (axial, z-direction)
    case 2
        t0 = 2; Tr = 3;
        Wmax = [0;0;4];
        if t >= t0 && t <= t0+Tr
            W_I = Wmax*(t-t0)/Tr;
        elseif t > t0+Tr
            W_I = Wmax;
        end

    %% 3. 1-Cosine gust (axial)
    case 3
        t0 = 3; Tg = 3;
        Wmax = [0;0;3];
        if t >= t0 && t <= t0+Tg
            W_I = 0.5*Wmax*(1-cos(pi*(t-t0)/Tg));
        elseif t > t0+Tg
            W_I = Wmax;
        end

    %% 4. Ramp with randomized direction in x-z plane
    case 4
        t0 = 2; Tr = 3;
        theta = pi/6; % fixed direction (avoid time-varying randomness)
        dir = [sin(theta);0;cos(theta)];
        Wmax = 4*dir;
        if t >= t0 && t <= t0+Tr
            W_I = Wmax*(t-t0)/Tr;
        elseif t > t0+Tr
            W_I = Wmax;
        end

    %% 5. Overlay ramp + cosine (random directions)
    case 5
        theta1 = pi/4;
        theta2 = -pi/6;

        dir1 = [sin(theta1);0;cos(theta1)];
        dir2 = [sin(theta2);0;cos(theta2)];

        % ramp
        t0=2; Tr=3;
        W1=[0;0;0];
        if t>=t0 && t<=t0+Tr
            W1=3*dir1*(t-t0)/Tr;
        elseif t>t0+Tr
            W1=3*dir1;
        end

        % cosine
        t1=6; Tg=2;
        W2=[0;0;0];
        if t>=t1 && t<=t1+Tg
            W2=2*dir2*0.5*(1-cos(pi*(t-t1)/Tg));
        elseif t>t1+Tg
            W2=2*dir2;
        end

        W_I=W1+W2;

    %% 6. Noisy overlay
    case 6
        W_I = gust_generator_overlay(t);
        noise_std = 0.1;
        W_I = W_I + noise_std*randn(3,1);

end

end
