clear all
%!Initialise constants
Pi = 3.141592653589793238462643;  %! In Matlab, can also simply use pi
g = 9.81;                         %! Acceleration due to gravity (m s^-2)
c_pa = 1005.0;                    %! Specific heat capacity of dry air (J kg^-1 K^-1)
Rho_w = 1000.0;                   %! Density of liquid water (Kg m^-3)
Rho_a = 1.225;			          %! Density of air (Kg m^-3)
Eps = 0.622;                      %! Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;                       %! Latent heat of vapourisation (J Kg^-1)
Ra = 287.0;                       %! Gas constant of dry air (J kg^-1 K^-1)
Rv = 462.0;                       %! Gas constant of water vapour (J kg^-1 K^-1)
k = 0.024;                        %! Thermal Conductivity of Air (J m^-1 s^-1 K^-1)
Kv = 2.21e-5 ;                    %! Diffusivity of Water Vapour (m^2 s^-1)
%es;		 		              %! Saturation vapor pressure (provided by function svp.m)

% Initial conditions

s = 0.003;        % constant supersaturation percentage
T = 282;          % constant temperature in Kelvin
r(1) = 1e-6;      % initial radius in metres

% Constant A3 from Davenish et al.
A3 = ((Lv^2 * Rho_w)/(k*Rv*T^2) + (Rho_w * Rv * T)/(Kv*svp(T)))^(-1);

r_star = 0.5e-3;
t = (r_star^2 - 2*c)/(2*A3*s) %! Time taken to grow (s)
