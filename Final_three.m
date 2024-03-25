function [R, S, T, P] = three(r, s, T, P)

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
w = 0.3;                          % Constant vertical velocity (m/s)
c = 1e+8;                         % CCN conc. at s
N = 100e6;

% Constants
Q1 = (1/T)*(((Eps*Lv*g)/(Ra*c_pa*T))-(g/Ra));
Q2 = Rho_a*(((Ra*T)/(Eps*svp(T)))+((Eps*(Lv^2))/(P*T*c_pa)));

% Calculate water vapor mixing ratio
qv = (Ra/Rv) * (svp(T) / P); 

% Calculate coefficients A1, A2, A3
A1 = (g / (Ra * T)) * (((Lv * Ra) / (c_pa * Rv * T)) - 1); 
A2 = ((Lv^2) / (c_pa * Rv * (T^2))) + (1 / qv);           
A3 = ((((Lv^2) * Rho_w) / (k * Rv * T^2)) + ((Rho_w * Rv * T) / (Kv * svp(T))))^-1;

P = -(P*g*w)/(Ra*T);
ql = ((4 * Pi * Rho_w * N) / Rho_a) * (r^2) * (A3 * (s / r));  % Calculate liquid water mixing ratio
T = (-(g*w)/(c_pa)) + (Lv/c_pa)*ql;

R = A3* (s/r);
S = (A1 * w) - (A2 * ql);
end
