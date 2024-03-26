clear all;

% Constants
g = 9.81;               % Acceleration due to gravity (m/s^2)
c_pa = 1005.0;          % Specific heat capacity of dry air (J/kg/K) 
Rho_w = 1000.0;         % Density of liquid water (kg/m^3)
Rho_a = 1.225;          % Density of air (kg/m^3)
Eps = 0.62;             % Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;             % Latent heat of vaporisation (J/kg) 
Ra = 287.0;             % Gas constant of dry air (J/kg/K)
Rv = 462.0;             % Gas constant of water vapour (J/kg/K)
k = 0.024;              % Thermal Conductivity of Air (J/m/s/K)
Kv = 2.21e-5;           % Diffusivity of Water Vapour (m^2/s)

%Initial Values
cloud_base = 250;                                   % Cloud depth in meters
initial_R_values = linspace(1e-6, 1e-5, 100);       % Initial droplet radius values (m)
initial_N_values = linspace(500e5,500e6 , 100);     % Initial droplet number concentration values (cm^-3)
initial_w_values = linspace(0.1, 2, 100);           % Define the values of initial vertical velocity
dt = 0.01;                                          % Time step (s)
z = linspace(0, cloud_base, cloud_base/dt);         % Vector of heights from cloud base
n_steps = length(z) - 1;

% Arrays to store results
max_R_R = zeros(size(initial_R_values));            % Max droplet radius for varying initial radius
max_S_R = zeros(size(initial_R_values));            % Max supersaturation for varying initial radius
max_R_N = zeros(size(initial_N_values));            % Max droplet radius for varying initial N
max_S_N = zeros(size(initial_N_values));            % Max supersaturation for varying initial N
max_R_w = zeros(size(initial_w_values));            % Max droplet radius for varying initial vertical velocity
max_S_w = zeros(size(initial_w_values));            % Max supersaturation for varying initial vertical velocity
P = zeros(1, length(n_steps));
T = zeros(1, length(n_steps));
T(1) = 282;                                         % Initial temperature
P(1) = 100000;                                       % Initial pressure

%time stepping for varying initial radius
for r_index = 1:length(initial_R_values)
    N = 100e6;
    w = 0.3;                                        % Vertical velocity (m/s)
    r(1) = initial_R_values(r_index);               % Set initial droplet radius
    s(1) = 0.003;                                   % Initial supersaturation ratio
    P(1) = 100000;
    T(1) = 282;
    for n = 1:n_steps
        [k1_R, k1_S, k1_T, k1_P] = four(r(n), s(n), T(n), P(n), N, w);
        [k2_R, k2_S, k2_T, k2_P] = four(r(n) + (dt/2)*k1_R, s(n) + (dt/2)*k1_S, T(n) + (dt/2)*k1_T, P(n) + (dt/2)*k1_P, N, w);
        [k3_R, k3_S, k3_T, k3_P] = four(r(n) + (dt/2)*k2_R, s(n) + (dt/2)*k2_S, T(n) + (dt/2)*k2_T, P(n) + (dt/2)*k2_P, N, w);
        [k4_R, k4_S, k4_T, k4_P] = four(r(n) + dt*k3_R, s(n) + dt*k3_S, T(n) + dt*k3_T, P(n) + dt*k3_P, N, w);

        T(n+1) = T(n) + (dt/6) * (k1_T + 2*k2_T + 2*k3_T + k4_T);
        P(n+1) = P(n) + (dt/6) * (k1_P + 2*k2_P + 2*k3_P + k4_P);

        % Update variables using RK4
        r(n+1) = r(n) + (dt/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
        s(n+1) = s(n) + (dt/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
    end
    max_R_R(r_index) = max(r);                      % Store max droplet radius
    max_S_R(r_index) = max(s);                      % Store max supersaturation
end

%time stepping for varying initial droplet number concentration
for N_index = 1:length(initial_N_values)
    N = initial_N_values(N_index);                  % Set initial droplet number concentration
    s(1) = 0.003;                                   % Initial supersaturation ratio
    w = 0.3;                                        % Vertical velocity (m/s)
    r(1) = 1e-6; 
    P(1) = 100000;
    T(1) = 282;
    for n = 1:n_steps
        [k1_R, k1_S, k1_T, k1_P] = four(r(n), s(n), T(n), P(n), N, w);
        [k2_R, k2_S, k2_T, k2_P] = four(r(n) + (dt/2)*k1_R, s(n) + (dt/2)*k1_S, T(n) + (dt/2)*k1_T, P(n) + (dt/2)*k1_P, N, w);
        [k3_R, k3_S, k3_T, k3_P] = four(r(n) + (dt/2)*k2_R, s(n) + (dt/2)*k2_S, T(n) + (dt/2)*k2_T, P(n) + (dt/2)*k2_P, N, w);
        [k4_R, k4_S, k4_T, k4_P] = four(r(n) + dt*k3_R, s(n) + dt*k3_S, T(n) + dt*k3_T, P(n) + dt*k3_P, N, w);

        T(n+1) = T(n) + (dt/6) * (k1_T + 2*k2_T + 2*k3_T + k4_T);
        P(n+1) = P(n) + (dt/6) * (k1_P + 2*k2_P + 2*k3_P + k4_P);

        % Update variables using RK4
        r(n+1) = r(n) + (dt/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
        s(n+1) = s(n) + (dt/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
    end
    max_R_N(N_index) = max(r);                      % Store max droplet radius
    max_S_N(N_index) = max(s);                      % Store max supersaturation
end

%time stepping for varying initial velocity
for w_index = 1:length(initial_w_values)
    w = initial_w_values(w_index);                          % Set initial vertical velocity
    N = 100e6;
    r(1) = 1e-6;                                    % Initial droplet radius
    s(1) = 0.003;                                   % Initial supersaturation ratio
    P(1) = 100000;
    T(1) = 282;
    for n = 1:n_steps
        [k1_R, k1_S, k1_T, k1_P] = four(r(n), s(n), T(n), P(n), N, w);
        [k2_R, k2_S, k2_T, k2_P] = four(r(n) + (dt/2)*k1_R, s(n) + (dt/2)*k1_S, T(n) + (dt/2)*k1_T, P(n) + (dt/2)*k1_P, N, w);
        [k3_R, k3_S, k3_T, k3_P] = four(r(n) + (dt/2)*k2_R, s(n) + (dt/2)*k2_S, T(n) + (dt/2)*k2_T, P(n) + (dt/2)*k2_P, N, w);
        [k4_R, k4_S, k4_T, k4_P] = four(r(n) + dt*k3_R, s(n) + dt*k3_S, T(n) + dt*k3_T, P(n) + dt*k3_P, N, w);

        T(n+1) = T(n) + (dt/6) * (k1_T + 2*k2_T + 2*k3_T + k4_T);
        P(n+1) = P(n) + (dt/6) * (k1_P + 2*k2_P + 2*k3_P + k4_P);

        % Update variables using RK4
        r(n+1) = r(n) + (dt/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
        s(n+1) = s(n) + (dt/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
    end
    max_R_w(w_index) = max(r);                      % Store max droplet radius
    max_S_w(w_index) = max(s);                      % Store max supersaturation
end

% Plotting
% Initial Vertical Velocity vs Max Droplet Radius
figure;
subplot(2, 1, 1);
plot(initial_w_values, max_R_w);
xlabel('Initial Vertical Velocity (m/s)');
ylabel('Max Droplet Radius (m)');
title('Max Droplet Radius vs Initial Vertical Velocity');

% Initial Vertical Velocity vs Max Supersaturation
subplot(2, 1, 2);
plot(initial_w_values, max_S_w);
xlabel('Initial Vertical Velocity (m/s)');
ylabel('Max Supersaturation');
title('Max Supersaturation vs Initial Vertical Velocity');

% Initial Droplet Radius vs Max Droplet Radius
figure;
subplot(2, 1, 1);
plot(initial_R_values, max_R_R);
xlabel('Initial Droplet Radius (m)');
ylabel('Max Droplet Radius (m)');
title('Max Droplet Radius vs Initial Droplet Radius');

% Initial Droplet Radius vs Max Supersaturation
subplot(2, 1, 2);
plot(initial_R_values, max_S_R);
xlabel('Initial Droplet Radius (m)');
ylabel('Max Supersaturation');
title('Max Supersaturation vs Initial Droplet Radius');

% Initial Droplet Number Concentration vs Max Droplet Radius
figure;
subplot(2, 1, 1);
plot(initial_N_values, max_R_N);
xlabel('Initial Droplet Number Concentration (cm^{-3})');
ylabel('Max Droplet Radius (m)');
title('Max Droplet Radius vs Initial Droplet Number Concentration');

% Initial Droplet Number Concentration vs Max Supersaturation
subplot(2, 1, 2);
plot(initial_N_values, max_S_N);
xlabel('Initial Droplet Number Concentration (cm^{-3})');
ylabel('Max Supersaturation');
title('Max Supersaturation vs Initial Droplet Number Concentration')
