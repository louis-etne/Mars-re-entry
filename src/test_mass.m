addpath('Data');
load('constants.mat')

peak_heatflux = 10.7; %kwh/m^2 : obtained with simulink simulation (formula slide 11 course 2)
r_n = 1.59; % obtained with Catia Model and Nasa official data


% Method 1 :
fun = @(x) (0.24 * (peak_heatflux.*cosd(x).^1.2) + 0.29 * sqrt(peak_heatflux.*cosd(x).^1.2) + 11.3) .* sind(x); % formula mass ablative shield

delta_max = 40;
m_tot = r_n^2 * 2 * pi * integral(fun, 0, delta_max); % obtain with Catia model and Nasa official documentation

% Method 2 : q_stag appy on every point of the capsule

heatshield_specific_mass = 0.24*peak_heatflux+0.29*sqrt(peak_heatflux)+11.3;
surface = 23.226; % Total surface calculated under catia

Mass_heatshield2 = heatshield_specific_mass * surface;

% Method 3 : plateau
mass_tot = 0;
zones = [0, 0.1, 0.3, 0.5, 0.8, 1];

for i = 1 : length(zones) - 1
    phi = delta_max * zones(i);
    mass_needed = 0.24 * (peak_heatflux.*cosd(phi).^1.2) + 0.29 * sqrt(peak_heatflux.*cosd(phi).^1.2) + 11.3;
    mass_tot = mass_tot + mass_needed * Vehicle.S * (zones(i + 1) - zones(i));
end

mass_tot = mass_tot + mass_needed * (surface - Vehicle.S);

