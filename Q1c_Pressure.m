% Constants
clear all
Pi = 3.141592653589793238462643;   % Pi
g = 9.81;                          % Acceleration due to gravity (m/s^2)
c_pa = 1005.0;                     % Specific heat capacity of dry air (J/kg/K)
Rho_w = 1000.0;                    % Density of liquid water (kg/m^3)
Rho_a = 1.225;                     % Density of air (kg/m^3)
Eps = 0.622;                       % Ratio of molecular masses of water vapor to dry air
Lv = 2.5e6;                        % Latent heat of vaporization (J/kg)
Ra = 287.0;                        % Gas constant for dry air (J/kg/K)
Rv = 462.0;                        % Gas constant for water vapor (J/kg/K)
k = 0.024;                         % Thermal conductivity of air (J/m/s/K)
Kv = 2.21e-5;                      % Diffusivity of water vapor (m^2/s)

% Initial radii in meters (converted from micrometers)
radii_initial = [0.01e-6, 0.1e-6, 1e-6, 10e-6, 100e-6]; 

s = 0.003;             % Constant supersaturation percentage
T = 282;               % Constant temperature in Kelvin

% Set timestep and number of steps
dt = 0.1;              % Timestep in seconds
n_steps = 45*60/dt;    % Number of steps for 45 minutes

% Time array
t = 0:dt:45*60;         % Time up to 2700 seconds (45 minutes)

% Initialize matrix for storing radius data
radius_matrix = zeros(length(t), length(radii_initial));

% Pre-calculate A3 using the constant temperature
A3 = ((Lv^2 * Rho_w)/(k*Rv*T^2) + (Rho_w * Rv * T)/(Kv*svp(T)))^(-1);

% Loop through each initial radius
for idx = 1:length(radii_initial)
    r_initial = radii_initial(idx);
    c = 0.5*r_initial^2;
    rstar = sqrt(2*A3*s*t + 2*c);
    radius_matrix(:, idx) = rstar * 10^6; % Convert to microns for plotting
end

% Plotting
figure;
hold on;
colors = ['b', 'g', 'r', 'c', 'm']; % Colour for each line
for idx = 1:length(radii_initial) %Plotting graphs
    plot(t, radius_matrix(:, idx), 'DisplayName', [num2str(radii_initial(idx)*10^6), '\mum']);
end
hold off;
legend('show');
xlabel('Time (s)');
ylabel('Droplet Radius (\mu m)');
title('Growth of Droplets Over Time at 282K for Various Initial Radii');
%ylim([0,20]) To zoom in 
%xlim([0,500])