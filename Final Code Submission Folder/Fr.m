function f = Fr(r)

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
T = 282;
s=0.003;

A3 = ((Lv^2*Rho_w)/(k*Rv*T^2) + (Rho_w*Rv*T)/(Kv*svp(T)))^-1;

f = A3 * (s / r);
end