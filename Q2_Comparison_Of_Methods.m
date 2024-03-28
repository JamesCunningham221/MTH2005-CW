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

% Calculate water vapor mixing ratio
qv = (Ra/Rv) * (svp(T) / P); 

% Calculate coefficients A1, A2, A3
A1 = (g / (Ra * T)) * (((Lv * Ra) / (c_pa * Rv * T)) - 1); 
A2 = ((Lv^2) / (c_pa * Rv * (T^2))) + (1 / qv);           
A3 = ((((Lv^2) * Rho_w) / (k * Rv * T^2)) + ((Rho_w * Rv * T) / (Kv * svp(T))))^-1;

 N = 100e6  ; % Droplet number density

% Time step (s)
dt = 12; 
cloud_base = 400;
% Vector of heights from cloud base
z = linspace(0,cloud_base, cloud_base/dt);
n_steps = length(z) - 1;

% Arrays to store results
s = zeros(1, n_steps+1);             % Supersaturation
r = zeros(1, n_steps+1);             % Droplet radius
s_percent = zeros(1, n_steps+1);     % Supersaturation

% Initial Conditions
r(1) = 1e-6;         % Initial droplet radius (m)
s(1) = 0;            % Initial supersaturation ratio


% Runge-Kutta integration
for n = 1:n_steps
    [k1_R, k1_S] = two(r(n), s(n));
    [k2_R, k2_S] = two(r(n) + (dt/2)*k1_R, s(n) + (dt/2)*k1_S);
    [k3_R, k3_S] = two(r(n) + (dt/2)*k2_R, s(n) + (dt/2)*k2_S);
    [k4_R, k4_S] = two(r(n) + dt*k3_R, s(n) + dt*k3_S);

    % Update variables using RK4
    r(n+1) = r(n) + (dt/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
    s(n+1) = s(n) + (dt/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);

    % Calculate supersaturation percentage
    s_percent(n+1) = s(n+1)* 100;
end

% Plot Supersaturation against altitude RK4
subplot(1, 2, 1);
plot(s_percent, z, 'b');
xlabel('Supersaturation');
ylabel('Altitude above cloud base (m)');
title('Supersaturation Percentage Evolution');
hold on

% Plot Droplet Radius against altitude RK4
subplot(1, 2, 2);
plot(r, z, 'b');
xlabel('Droplet Radius (m)');
ylabel('Altitude (m)');
title('Droplet Radius Evolution');

hold on

% Forward Euler time stepping
for i = 1:n_steps

    % Work out derivatives of r and q_w with respect to t
    dr_dt = A3 * s(i)/r(i);
    dql_dt = (4 * pi * Rho_w * N / Rho_a) * r(i)^2 * dr_dt;

    % Work out r and s
    r(i+1) = r(i) + dt * dr_dt;
    s(i+1) = s(i) + dt * (A1 * w - A2 * dql_dt);
    % Calculate supersaturation percentage
    s_percent(i+1) = s(i+1)* 100;
end

% Plot Forward Euler
% figure;
% 

% Plot Supersaturation against altitude FE
subplot(1, 2, 1);
plot(s_percent, z, 'r');
xlabel('Supersaturation (%)');
ylabel('Altitude above cloud base (m)');
title('Supersaturation Percentage Evolution');
legend('RK4', 'Forward-Euler')

% Plot Droplet Radius against altitude FE
subplot(1, 2, 2);
plot(r, z, 'r');
xlabel('Droplet Radius (m)');
ylabel('Altitude (m)');
title('Droplet Radius Evolution');
legend('RK4','Forward-Euler')




hold off