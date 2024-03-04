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

r(1) = 1e-10;                        %really small number
s(1) = 1e-10;                        %supersaturation ratio
N = 5;                        %droplet number density
w = 1000;                          %initial altitude
qv = 0.0713 ;                      %water vapour mixing ratio
P = 10000;                       %Constant Pressure in Pascals

%EquationsT
%qv = 0.622*svp(T)/P; %water vaping mixing ratio
A1 = (g/(Ra*T))*((Lv*Ra)/(c_pa*Rv*T)-1);
A2 = ((Lv^2)/(c_pa*Rv*(T^2)))+(1/qv);
A3 = (((Lv^2)*Rho_w)/(k*Rv*T^2)) + ((Rho_w*Rv*T)/(Kv*svp(T)))^-1;
%s = (e/svp(T))-1 %supersaturation ratio
%dr = A3*(s/r); %droplet radius equation
%dql = ((4*pi*Rho_w*N)/(Rho_a))*(r^2)*(dr); %liquid water mixing ratio equation
%dz = w; %altitude equation
%ds = A1*w -A2*(dql); %supersaturation equation

% Time parameters
delt = 0.01;            % Time step (s)
t_end = 60;      % Total simulation time (s)
nsteps = t_end / delt; % Number of time steps

% Arrays to store results
altitude = zeros(nsteps, 1);
supersaturation = zeros(nsteps, 1);
droplet_radius = zeros(nsteps, 1);
%liquid_water_mixing_ratio = zeros(nsteps, 1);

subplot(1, 2, 1);


% Forward Euler time stepping for radius
for n = 1:nsteps
    % Calculate water vapor mixing ratio
    qv = Eps * svp(T) / P;



    % Calculate droplet radius
    dr = A3 * (s(n) / r(n));

    % Calculate liquid water mixing ratio
    dql = ((4 * Pi * Rho_w * N) / Rho_a) * (r(n)^2) * dr;

    % Update altitude
    dz = w(n);

    % Update supersaturation
    ds = A1 * w(n) - A2 * dql;

    % Update droplet radius, liquid water mixing ratio, and altitude
    r(n+1) = r(n) + dr * delt;
    %N = N + dql * delt;
    w(n+1) = w(n) + dz * delt ;
    s(n+1) = s(n) + ds * delt;

    % Store results
    altitude(n) = w(n);
    supersaturation(n) = s(n);
    droplet_radius(n) = r(n);
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
plot(altitude, supersaturation);
xlabel('Altitude (m)');
ylabel('Supersaturation');
title('Supersaturation Evolution');
% Zoom in on a specific region of the plot
%xlim([0, 1000]); % Adjust the limits according to desired region
%ylim([0, 0.5]); % Adjust the limits according to desired region

subplot(1, 2, 2);
plot(altitude, droplet_radius);
xlabel('Altitude (m)');
ylabel('Droplet Radius (m)');
title('Droplet Radius Evolution');
