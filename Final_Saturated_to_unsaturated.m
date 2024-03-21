clearvars -except Pi g c_pa Rho_w Rho_a Eps Lv Ra Rv k Kv
r0 = 7e-6;                  % Initial radius of droplet in microns
s = 75/100 -1;              % Supersaturation in percentage
T = 273:0.25:373;           % Possible temperatures of water (0 degrees Celsius to 100 degrees Celsius) in vector
dt = 0.1;                   % Timestep
d = 2.45e-11;               % Constant value when t = 0


for i = 1:length(T)
    A3 = ((Lv^2 * Rho_w)/(k*Rv*T(i)^2) + (Rho_w * Rv * T(i))/(Kv*svp(T(i))))^(-1);
    t(i) = -(r0^2 - d) / (2*A3*s);
end

% Plot time taken vs temperature
plot(T, t)
xlabel('Temperature [Kelvin , K]') 
ylabel('Time [seconds , s]')
title({'Time Taken for an Isolated Cloud Droplet to Evaporate Completely'}) 
subtitle('by Water Vapour Diffusion')
grid on
