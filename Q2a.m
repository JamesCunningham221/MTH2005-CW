clear

Pi = 3.141592653589793238462643;  %! In Matlab, can also simply use pi
g = 9.81;                         %! Acceleration due to gravity (m s^-2)
c_pa = 1005.0;                    %! Specific heat capacity of dry air (J kg^-1 K^-1)
Rho_w = 1000.0;                   %! Density of liquid water (Kg m^-3)
Rho_a = 1.225; 			         %! Density of air (Kg m^-3)
Eps = 0.622;                      %! Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;                       %! Latent heat of vapourisation (J Kg^-1)
Ra = 287.0;                       %! Gas constant of dry air (J kg^-1 K^-1)
Rv = 462.0;                       %! Gas constant of water vapour (J kg^-1 K^-1)
k = 0.024;                        %! Thermal Conductivity of Air (J m^-1 s^-1 K^-1)
Kv = 2.21e-5;                     %! Diffusivity of Water Vapour (m^2 s^-1)
T = 282;                         %Constant Temperature in Kelvin

% Time parameters
dt = 1;            % Time step (s)
t_end = 90*60;      % Total simulation time (s)
nsteps = t_end / dt; % Number of time steps
t = 0:dt:nsteps;

% Arrays to store results
z = zeros(1, length(t)); %altitude
s = zeros(1, length(t)); %supersaturation
r = zeros(1, length(t)); %droplet radius
%liquid_water_mixing_ratio = zeros(nsteps, 1);

r(1) = 1e-6;                        %really small number
s(1) = 1.5;                        %supersaturation ratio
N = 100e-6;                        %droplet number density
z(1) = 1000;                          %initial altitude 
%qv = 0.0713 ;                      %water vapour mixing ratio
P = 100000;                       %Constant Pressure in Pascals
w = 1;
 % Calculate water vapor mixing ratio
qv = (Ra/Rv) * svp(T) / P;
%EquationsT
%qv = 0.622*svp(T)/P; %water vaping mixing ratio
A1 = (g/(Ra*T))*(((Lv*Ra)/(c_pa*Rv*T))-1);
A2 = ((Lv^2)/(c_pa*Rv*(T^2)))+(1/qv);
A3 = ((((Lv^2)*Rho_w)/(k*Rv*T^2)) + ((Rho_w*Rv*T)/(Kv*svp(T))))^-1;
%s = (e/svp(T))-1 %supersaturation ratio
%dr = A3*(s/r); %droplet radius equation
%dql = ((4*pi*Rho_w*N)/(Rho_a))*(r^2)*(dr); %liquid water mixing ratio equation
%dz = w; %altitude equation
%ds = A1*w -A2*(dql); %supersaturation equation

subplot(1, 2, 1)  

% Forward Euler time stepping for radius
for n = 1:nsteps
    % Calculate droplet radius
    r(n+1) = r(n) + dt*(A3 * (s(n)./r(n)));

    % Calculate liquid water mixing ratio
    ql = ((4 * Pi * Rho_w * N) / Rho_a) * (r(n).^2) * r(n+1);

    % Update altitude
    z(n+1) = z(n) + dt*w;

    % Update supersaturation
    s(n+1) = s(n) + dt*(A1 * w - A2 * ql);
end


% % Forward Euler time stepping for supersaturation
% for n = 1:nsteps
%     % Calculate water vapor mixing ratio
%     qv = Eps * svp(T) / P;
% 
%     % Calculate coefficients for equations
%     A1 = (g / (Ra * T)) * ((Lv * Ra) / (c_pa * Rv * T) - 1);
%     A2 = ((Lv^2) / (c_pa * Rv * (T^2))) + (1 / qv);
%     A3 = (((Lv^2) * Rho_w) / (k * Rv * (T^2))) + ((Rho_w * Rv * T) / (Kv * svp(T)));
% 
%     % Calculate droplet radius
%     dr = A3 * (s(n) / r);
% 
%     % Calculate liquid water mixing ratio
%     dql = ((4 * Pi * Rho_w * N) / Rho_a) * (r^2) * dr;
% 
%     % Update altitude
%     dz = w(n);
% 
%     % Update supersaturation/
%     ds = A1 * w - A2 * dql;
% 
%     % Update droplet radius, liquid water mixing ratio, and altitude
%     %r = r + dr * delt;
%     %N = N + dql * delt;
%     w(n+1) = w(n) + dz * delt ;
%     s(n+1) = s(n) + ds * delt;
% 
%     % Store results
%     altitude(n) = w(n);
%     supersaturation(n) = s(n);
%     %droplet_radius(t) = r;
% end

% Plotting
figure;
subplot(1, 2, 1);
plot(z, s);
xlabel('Altitude (m)');
ylabel('Supersaturation');
title('Supersaturation Evolution');
% Zoom in on a specific region of the plot
%xlim([0, 1000]); % Adjust the limits according to desired region
%ylim([0, 0.5]); % Adjust the limits according to desired region

subplot(1, 2, 2);
plot(z, r);
xlabel('Altitude (m)');
ylabel('Droplet Radius (m)');
title('Droplet Radius Evolution');