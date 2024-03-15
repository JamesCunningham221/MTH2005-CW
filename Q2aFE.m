% Clearing any existing variables
clear

% Constants
Pi = 3.141592653589793238462643;  % Pi
g = 9.81;                         % Acceleration due to gravity (m/s^2)
c_pa = 1005.0;                    % Specific heat capacity of dry air (J/kg/K)
Rho_w = 1000.0;                   % Density of liquid water (kg/m^3)
Rho_a = 1.225; 			          % Density of air (kg/m^3)
Eps = 0.622;                      % Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;                       % Latent heat of vaporisation (J/kg)
Ra = 287.0;                       % Gas constant of dry air (J/kg/K)
Rv = 462.0;                       % Gas constant of water vapour (J/kg/K)
k = 0.024;                        % Thermal conductivity of air (J/m/s/K)
Kv = 2.21e-5;                     % Diffusivity of water vapour (m^2/s)
T = 282;                          % Constant temperature in Kelvin
P = 100000;                       % Constant pressure (Pa)
w = 0.3;                          % Constant vertical velocity (m/s)

% Constants
Q1 = (1/T)*(((Eps*Lv*g)/(Ra*c_pa*T))-(g/Ra));
Q2 = Rho_a*(((Ra*T)/(Eps*svp(T)))+((Eps*(Lv^2))/(P*T*c_pa)));
%Q1 = 5.82e-4;      % Production term coefficient (m^-1) at 273K
%Q2 = 423.6013;     % Condensation term coefficient (kg/kg) at 273K and 1000hPa

% Time parameters
dt = 0.1;            % Time step (s)
t_end = 0.5e5;     % Total simulation time (s)
nsteps = t_end / dt; % Number of time steps

% Arrays to store results
z = zeros(1, nsteps+1);            % Altitude
s = zeros(1, nsteps+1);            % Supersaturation
r = zeros(1, nsteps+1);            % Droplet radius
s_percentage = zeros(1, nsteps+1); % Supersaturation percentage

% Initial conditions
r(1) = 0.1e-3;    % Initial droplet radius (m)
s(1) = 0;     % Initial supersaturation ratio
z(1) = 0;    % Initial altitude (m)


% Calculate water vapor mixing ratio
qv = (Ra/Rv) * (svp(T) / P); 

% Calculate coefficients A1, A2, A3
A1 = (g / (Ra * T)) * (((Lv * Ra) / (c_pa * Rv * T)) - 1); 
A2 = ((Lv^2) / (c_pa * Rv * (T^2))) + (1 / qv);           
A3 = ((((Lv^2) * Rho_w) / (k * Rv * T^2)) + ((Rho_w * Rv * T) / (Kv * svp(T))))^-1;

% Forward Euler time stepping
for n = 1:nsteps

    N = c_pa*(s(n)*100)^k;     % Droplet number density

    % Calculate droplet radius
    r(n+1) = r(n) + dt * (A3 * (s(n) / r(n)));
    
    % Calculate liquid water mixing ratio
    ql = ((4 * Pi * Rho_w * N) / Rho_a) * (r(n)^2) * (A3 * (s(n) / r(n)));
    
    % Update altitude
    z(n+1) = z(n) + (dt * w);
    
    % Update supersaturation
    s(n+1) = s(n) + dt * ((A1 * w) - (A2 * ql));

    % Calculate supersaturation percentage
    s_percentage(n+1) = (s(n+1) - 1) * 100;
end

% Plot Supersaturation against altitude
subplot(1, 2, 1);
plot(s, z);
xlabel('Supersaturation (%)');
ylabel('Altitude above cloud base (m)');
title('Supersaturation Percentage Evolution');

% Plot Droplet Radius against altitude
subplot(1, 2, 2);
plot(r, z);
xlabel('Droplet Radius (m)');
ylabel('Altitude (m)');
title('Droplet Radius Evolution');
