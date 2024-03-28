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

% Constants for simulation
s = 0.003;                      % Supersaturation (0.30%)
T = 282;                        % Temperature (K)
dt = 0.1;                        % Time step (seconds)
num_steps = 45*60;              % Total number of steps (45 minutes * 60 seconds)

% Initialize arrays
t = 0:dt:num_steps;
r = zeros(1, length(t));        % Droplet radius
r_FE = zeros(1, length(t));        % Droplet radius

t(1) = 0;
r(1) = 1e-6;
r_FE(1) = 1e-6;
A3 = ((Lv^2*Rho_w)/(k*Rv*T^2) + (Rho_w*Rv*T)/(Kv*svp(T)))^-1;

for i = 1:length(t)
   k1 = Fr(r(i));
   k2 = Fr(r(i) + 0.5*dt*k1);
   k3 = Fr(r(i) + 0.5*dt*k2);
   k4 = Fr(r(i) + dt*k3);
   r(i+1) = r(i) + 1/6*dt*(k1 + 2*k2 + 2*k3 + k4);
end

for i = 1:length(t)
    % Calculate rate of change of droplet radius    
    % Update droplet radius
   r_FE(i+1) = r_FE(i) + dt*(A3 * (s / r_FE(i)));
end

% Calculate the true solution
c = 0.5*r(1)^2;
rstar = (2*A3*s*t + 2*c).^(1/2);

for n = 1:length(t)
    fe_err(n) = abs((r_FE(n) - rstar(n)));
    rk_err(n) = abs((r(n) - rstar(n)));
end


figure(1);
loglog(t, fe_err, 'b', 'LineWidth', 2)
hold on
loglog(t, rk_err, 'r', 'LineWidth', 2)
xlabel('Time (s)', 'FontSize', 15);
ylabel('Absolute Error', 'FontSize', 15);
title('Absolute errors of numerical schemes over time', 'FontSize', 15)
legend('Forward Euler', 'RK4', 'FontSize', 15)