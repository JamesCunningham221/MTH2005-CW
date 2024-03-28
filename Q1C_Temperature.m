clear all
%constants
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
s = 0.003;

temperatures = [263, 273, 283, 293, 303]; % temperatures in Kelvin

% Set timestep and number of steps
dt = 0.1;                % timestep in seconds
n_steps = 45*60/dt;     % number of steps for 45 minutes

% Time array
t = 0:dt:2700;           % time up to 2700 seconds (45 minutes)

% Keep initial radius constant at 1 micron
r_initial = 1e-6; % initial radius in meters

% Initialize matrix for storing radius data for different temperatures
radius_matrix_temp = zeros(length(t), length(temperatures));

% Loop through each temperature
for idx = 1:length(temperatures)
    T_current = temperatures(idx);
    
    % Update A3 based on current temperature, including call to svp(T_current)
    A3_current = ((Lv^2 * Rho_w)/(k*Rv*T_current^2) + (Rho_w * Rv * T_current)/(Kv*svp(T_current)))^(-1);
    
    % Calculate growth for current temperature
    c = 0.5*r_initial^2;
    rstar = (2*A3_current*s*t + 2*c).^(1/2);
    radius_matrix_temp(:, idx) = rstar*10^6; % Convert to microns for plotting
end

% Plotting the results
figure; 
hold on;

% Iterate through each temperature for plotting
for idx = 1:length(temperatures) %Plots graphs
    plot(t, radius_matrix_temp(:, idx), 'DisplayName', [num2str(temperatures(idx)), 'K']);
end

hold off;
legend('show');
xlabel('Time (s)');
ylabel('Temperature K');
title('Growth of a 1 Micron Droplet Over Time at Various Temperatures');

