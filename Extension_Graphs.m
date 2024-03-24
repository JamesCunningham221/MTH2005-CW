clear all

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
cloud_base = 1000;       % Cloud depth in metres
w = 0.3;                % Vertical velocity
c = 1e+8;                % CCN conc. at s
N = 100e6;              %100 cm^3 density in SI

% Time step (s)
dt = 0.1; 

% Vector of heights from cloud base
z = linspace(0,cloud_base, cloud_base/dt);
n_steps = length(z) -1;

% Arrays to store results
s = zeros(1, n_steps+1);             % Supersaturation
r = zeros(1, n_steps+1);             % Droplet radius
s_percent = zeros(1, n_steps+1);     % Supersaturation

% Initial Conditions and Vector Initialization
initial_R = linspace(1e-7,1e-5, n_steps+1); % Initial radii vector
final_R = zeros(1, length(initial_R)); % Final radii after simulation

% Initial Conditions
r(1) = 1e-6;         % Initial droplet radius (m)
s(1) = 0;            % Initial supersaturation ratio
es = svp(T);
% Calculate water vapour mixing ratio
qv = (Ra/Rv) * (es/P);

% Calculate coefficients A1, A2, A3
A1 = (g / (Ra * T)) * (((Lv * Ra) / (c_pa * Rv * T)) - 1); 
A2 = ((Lv^2) / (c_pa * Rv * (T^2))) + (1 / qv);           
A3 = ((((Lv^2) * Rho_w) / (k * Rv * T^2)) + ((Rho_w * Rv * T) / (Kv * svp(T))))^-1;

% Forward Euler time stepping
for i = 1:length(initial_R)
    r(1) = initial_R(i); % Set initial droplet radius
    for n = 1:n_steps -1
        
        [k1_R, k1_S, k1_T, k1_P] = three(r(n), s(n), T(n), P(n));
        [k2_R, k2_S, k2_T, k2_P] = three(r(n) + (dt/2)*k1_R, s(n) + (dt/2)*k1_S, T(n) + (dt/2)*k1_T, P(n) + (dt/2)*k1_P);
        [k3_R, k3_S, k3_T, k3_P] = three(r(n) + (dt/2)*k2_R, s(n) + (dt/2)*k2_S, T(n) + (dt/2)*k2_T, P(n) + (dt/2)*k2_P);
        [k4_R, k4_S, k4_T, k4_P] = three(r(n) + dt*k3_R, s(n) + dt*k3_S, T(n) + dt*k3_T, P(n) + dt*k3_P);
        T(n+1) = T(n) + (dt/6) * (k1_T + 2*k2_T + 2*k3_T + k4_T);
        P(n+1) = P(n) + (dt/6) * (k1_P + 2*k2_P + 2*k3_P + k4_P);
        r(n+1) = r(n) + (dt/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
        s(n+1) = s(n) + (dt/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
        s_percent(n+1) = s(n) * 100;
    end
    final_R(i) = max(r); 
end

figure;
plot(initial_R, final_R);
xlabel('r 0');
ylabel('R');
title('R against r0');

