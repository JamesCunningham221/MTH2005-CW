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

clearvars -except Pi g c_pa Rho_w Rho_a Eps Lv Ra Rv k Kv
r0 = 7e-6;                  % Initial radius of droplet in microns
s = 75/100 -1;              % Supersaturation in percentage
T = 273:0.25:373;           % Possible temperatures of water (0 degrees Celsius to 100 degrees Celsius) in vector
dt = 0.1;                   % Timestep
d = 2.45e-11;               % Constant value when t = 0
rfinal = 0;

for i = 1:length(T) %Calc time it takes to evaporate
    A3 = ((Lv^2 * Rho_w)/(k * Rv * T(i)^2) + (Rho_w * Rv * T(i))/(Kv * svp(T(i))))^(-1);
    t(i) = (rfinal^2 - r0^2) / (2 * A3 *s);
end



% Plot time taken vs temperature
plot(T, t)
xlabel('Temperature [Kelvin , K]') 
ylabel('Time [seconds , s]')
title({'Time Taken for an Isolated Cloud Droplet to Evaporate Completely'}) 
subtitle('by Water Vapour Diffusion')
grid on
