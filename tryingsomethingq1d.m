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

% Time parameters
dt = 0.1;            % Time step (s)
t_end = 2*24*60*60;     % Total simulation time (s)
nsteps = t_end / dt; % Number of time steps

% Arrays to store results
z = zeros(1, nsteps+1);            % Altitude
s = zeros(1, nsteps+1);            % Supersaturation
r = zeros(1, nsteps+1);            % Droplet radius
s_percentage = zeros(1, nsteps+1); % Supersaturation percentage

% Initial conditions
r(1) = 1e-6;    % Initial droplet radius (m)
s(1) = 0.003;     % Initial supersaturation ratio
z(1) = 0;    % Initial altitude (m)
P(1) = 100000; % Constant pressure (Pa)
w = 0.3;          % Constant vertical velocity (m/s)
T(1) = 282;

% Calculate water vapor mixing ratio
qv = (Ra/Rv) * (svp(T) / P); 

% Forward Euler time stepping
for n = 1:nsteps

    % Calculate coefficients A1, A2, A3
A1 = (g / (Ra * T(n))) * (((Lv * Ra) / (c_pa * Rv * T(n))) - 1); 
A2 = ((Lv^2) / (c_pa * Rv * (T(n)^2))) + (1 / qv);           
A3 = ((((Lv^2) * Rho_w) / (k * Rv * T(n)^2)) + ((Rho_w * Rv * T(n)) / (Kv * svp(T(n)))))^-1;


    N = c_pa*(s(n)*100)^k;     % Droplet number density

    qv = (Ra/Rv) * (svp(T(n)) ./ P(n));

    % Calculate droplet radius
    drdt_FE = A3 * (s(n) / r(n));
    
    % Update droplet radius
    r(n+1) = r(n) + dt * drdt_FE;
    
    % Calculate liquid water mixing ratio
    ql = ((4 * Pi * Rho_w * N) / Rho_a) * (r(n)^2) * (A3 * (s(n) / r(n)));
    
    P(n+1) = P(n) + dt*((-g*P(n)*w)/(Ra*T(n)));

    T(n+1) = T(n)+ dt*(((-g*w)/c_pa)) + ((Lv/c_pa)*ql);
    
    % Update supersaturation
    s(n+1) = s(n) + dt * ((A1 * w) - (A2 * ql));
end

t = 0 : dt : t_end;
r(end)
plot(t, r)
ylabel('radius')
xlabel('time')