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

% Set timestep and number of steps to iterate
dt = 0.1;                % timestep
n_steps = 250*60;     % number of steps based on timestep
t = 0:dt:n_steps;           % time up to 2700 seconds (45 minutes)
r = 1e-6:(100e-6 - 1e-6)/n_steps:100e-6;   % initial radius in microns


% Constant A3 from Davenish et al.
A3 = ((Lv^2 * Rho_w)/(k*Rv*T^2) + (Rho_w * Rv * T)/(Kv*svp(T)))^(-1);

% Calculate the true solution
c = 0.5*r(1)^2;
rstar = (2*A3*s*t + 2*c).^(1/2);

% Plot true solution
plot(rstar*10^(6), t/60)

ylabel('Time [minutes]') 
xlabel('Droplet Radius [microns, \mum]')
title({'True Solution for the Growth of the';'Droplet Radius Over 45 Minutes'}) 
subtitle('for Temperature 282K and Initial Radius 1 [\mum]')
grid on

r_star = 100e-6;
g = (r_star^2 - 2*c)/(2*A3*s);

g %in second!