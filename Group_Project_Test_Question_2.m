clear

% Constants

g = 9.81;               %! Acceleration due to gravity (m s^-2)
c_pa = 1005.0;          %! Specific heat capacity of dry air (J kg^-1 K^-1) 
Rho_w = 1000.0;         %! Density of liquid water (Kg m^-3)
Rho_a = 1.225;          %! Density of air (Kg m^-3)
Eps = 0.62;             %! Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;             %! Latent heat of vapourisation (J Kg^-1) 
Ra = 287.0;             %! Gas constant of dry air (J kg^-1 K^-1)
Rv = 462.0;             %! Gas constant of water vapour (J kg^-1 K^-1)
k = 0.024;              %! Thermal Conductivity of Air (J m^-1 s^-1 K^-1)
Kv = 2.21e-5;           %! Diffusivity of Water Vapour (m^2 s^-1)

% Constants
T = 282;                % Constant temperature in Kelvin
P = 100000;             % Constant pressure in pascals N = 1e+8; % droplet number density
es = svp(T);            % Saturation vapour pressure 
cloud_base = 400;       % Cloud depth in metres
w = 0.3;                % Vertical velocity

% Time step (s)
dt = 0.1; 

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

% Calculate water vapour mixing ratio
qv = (Ra/Rv) * (es/P);

% Calculate coefficients A1, A2, A3
A1 = (g / (Ra * T)) * (((Lv * Ra) / (c_pa * Rv * T)) - 1); 
A2 = ((Lv^2) / (c_pa * Rv * (T^2))) + (1 / qv);           
A3 = ((((Lv^2) * Rho_w) / (k * Rv * T^2)) + ((Rho_w * Rv * T) / (Kv * svp(T))))^-1;

% Forward Euler time stepping
for i = 1:n_steps

    N = 1e+8;   % Droplet number density

    % Work out derivatives of r and q_w with respect to t
    dr_dt = A3 * s(i)/r(i);
    dql_dt = (4 * pi * Rho_w * N / Rho_a) * r(i)^2 * dr_dt;

    % Work out r and s
    r(i+1) = r(i) + dt * dr_dt;
    s(i+1) = s(i) + dt * (A1 * w - A2 * dql_dt);
    % Calculate supersaturation percentage
    s_percentage(i+1) = (s(i+1) - 1) * 100;
end

% Plot Forward Euler
figure;

% Plot Supersaturation against altitude
%subplot(1, 2, 1);
%plot(s, z);
%xlabel('Supersaturation');
%ylabel('Altitude above cloud base (m)');
%title('Supersaturation Evolution');

% Plot Droplet Radius against altitude
subplot(1, 2, 2);
plot(r, z);
xlabel('Droplet Radius (m)');
ylabel('Altitude (m)');
title('Droplet Radius Evolution');


% Plot Supersaturation against altitude
subplot(1, 2, 1);
plot(s, z);
xlabel('Supersaturation (%)');
ylabel('Altitude above cloud base (m)');
title('Supersaturation Percentage Evolution');




