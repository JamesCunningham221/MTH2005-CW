% Clearing any existing variables
clear

% Constants (unchanged)
Pi = 3.141592653589793238462643;  % Pi
g = 9.81;                         % Acceleration due to gravity (m/s^2)
c_pa = 1005.0;                    % Specific heat capacity of dry air (J/kg/K)
Rho_w = 1000.0;                   % Density of liquid water (kg/m^3)
Rho_a = 1.225; 			          % Density of air (kg/m^3)
Eps = 0.622;                      % Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;                       % Latent heat of vaporization (J/kg)
Ra = 287.0;                       % Gas constant of dry air (J/kg/K)
Rv = 462.0;                       % Gas constant of water vapour (J/kg/K)
k = 0.024;                        % Thermal conductivity of air (J/m/s/K)
Kv = 2.21e-5;                     % Diffusivity of water vapour (m^2/s)
T = 282;                          % Constant temperature in Kelvin
P = 100000;                       % Constant pressure (Pa)

% Define vertical velocity values
w_values = [0.1, 0.3, 0.5];

% Time step (s)
dt = 0.01; 

% Vector of heights from cloud base
cloud_base = 400;                 % Cloud base height (m)
z = linspace(0, cloud_base, cloud_base/dt);
n_steps = length(z) - 1;

% Arrays to store results for each velocity value
s_values = zeros(length(w_values), n_steps+1);  % Supersaturation
r_values = zeros(length(w_values), n_steps+1);  % Droplet radius
s_percent = zeros(1, n_steps+1);     % Supersaturation

% Runge-Kutta integration for each velocity value
for v = 1:length(w_values)
    w = w_values(v);
    
    % Initialize variables
    r = zeros(1, n_steps+1); % Droplet radius
    s = zeros(1, n_steps+1); % Supersaturation
    r(1) = 1e-6;             % Initial droplet radius (m)
    s(1) = 0;                % Initial supersaturation ratio
    
    % Runge-Kutta integration
    for n = 1:n_steps
        [k1_R, k1_S] = two(r(n), s(n), w);
        [k2_R, k2_S] = two(r(n) + (dt/2)*k1_R, s(n) + (dt/2)*k1_S, w);
        [k3_R, k3_S] = two(r(n) + (dt/2)*k2_R, s(n) + (dt/2)*k2_S, w);
        [k4_R, k4_S] = two(r(n) + dt*k3_R, s(n) + dt*k3_S, w);

        % Update variables using RK4
        r(n+1) = r(n) + (dt/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
        s(n+1) = s(n) + (dt/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
        s_percent(n+1) = s(n+1)* 100;
    end
    
    % Store results
    s_values(v, :) = s_percent;
    r_values(v, :) = r;
end

% Plot Supersaturation against altitude
subplot(1, 2, 1);
hold on;
for v = 1:length(w_values)
    plot(s_values(v, :), z, 'DisplayName', sprintf('w = %.1f m/s', w_values(v)));
end
xlabel('Supersaturation');
ylabel('Altitude above cloud base (m)');
title('Supersaturation Percentage Evolution');
legend('Location', 'best');

% Plot Droplet Radius against altitude
subplot(1, 2, 2);
hold on;
for v = 1:length(w_values)
    plot(r_values(v, :), z, 'DisplayName', sprintf('w = %.1f m/s', w_values(v)));
end
xlabel('Droplet Radius (m)');
ylabel('Altitude above cloud base(m)');
title('Droplet Radius Evolution');
legend('Location', 'best');