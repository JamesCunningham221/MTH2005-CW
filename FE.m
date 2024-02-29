clear all
% Program to solve growth of monodisperse cloud droplet population in an ascending air parcel

%--------------------------------------------------------------------------
% Initialise constants
Pi = pi;                        % In MATLAB, you can use 'pi'
g = 9.81;                       % Acceleration due to gravity (m/s^2)
c_pa = 1005.0;                  % Specific heat capacity of dry air (J/kg/K)
Rho_w = 1000.0;                 % Density of liquid water (kg/m^3)
Rho_a = 1.225;                  % Density of air (kg/m^3)
Eps = 0.622;                    % Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;                     % Latent heat of vaporization (J/kg)
Ra = 287.0;                     % Gas constant of dry air (J/kg/K)
Rv = 462.0;                     % Gas constant of water vapour (J/kg/K)
k = 0.024;                      % Thermal conductivity of air (J/m/s/K)
Kv = 2.21e-5;                   % Diffusivity of water vapour (m^2/s)
% Define initial vertical velocity and droplet number concentration as parameters
% Vertical velocity (m/s)
W = 0;
% Droplet Number Conc. (#/m^3)
NDropletDensity = 0;

% Constants for simulation
s = 0.003;                      % Supersaturation (0.30%)
T = 282;                        % Temperature (K)
dt =60;                        % Time step (seconds)
num_steps = 45*60;              % Total number of steps (45 minutes * 60 seconds)

% Initialize arrays
t = 0:dt:num_steps;
r = zeros(1, length(t));        % Droplet radius
t(1) = 0;
r(1) = 1e-6;
A3 = ((Lv^2*Rho_w)/(k*Rv*T^2) + (Rho_w*Rv*T)/(Kv*svp(T)))^-1;

% Main loop for simulation
for i = 1:length(t)-1
    % Calculate rate of change of droplet radius
    drdt_FE = A3 * (s / r(i));
    
    % Update droplet radius
    r(i+1) = r(i) + dt * drdt_FE;
end

%--------------------------------------------------------------------------
% Plot the results
plot(t, r, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Droplet Radius (m)');
title('Droplet Growth Over Time');
grid on;%---------------------------------------------------------------------------