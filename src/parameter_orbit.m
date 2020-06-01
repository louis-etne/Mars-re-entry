clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat', 'Mars', 'Orbit', 'Vehicle');
load('altitude.mat', 'alt', 'downrange');
%% Parameter of orbit around mars : 
% HYPOTHESIS : 
%   - circular orbit at 500 km altitude
% Supposition for the calculation : Mars Atmosphere begin at 120 km high
% speed departure orbit
V_departure_orbit = sqrt(2 * (-Mars.mu / (2 * Orbit.radius) + Mars.mu / Orbit.radius));

r_p = Mars.radius + 0 ;   % km - altitude of transfer orbit periapsis
a_i = Orbit.radius;      % km - semi-major axis of initial orbit
e_i = 0;                 % eccentricity of initial orbit

a_t = (Orbit.radius + r_p)/2;                        % km - sma of transfer orbit
e_t = (Orbit.radius - r_p) / (Orbit.radius + r_p);   % eccentricity of transfer orbit
Va_t = sqrt(Mars.mu * (2/Orbit.radius - 1/a_t));     % km/s velocity at apoapsis of transfer orbit
Vp_t = sqrt(Mars.mu * (2/r_p - 1/a_t));              % km/s velocity at periapsis of transfer orbit
DeltaV = Va_t - V_departure_orbit
h = r_p * Vp_t;
p = h^2 / Mars.mu;
r_entry = Mars.radius + 120;
theta_entry = acos(p/(r_entry * e_t) - 1/e_t);
V_entry = sqrt(Mars.mu * (2/r_entry - 1/a_t));

InitialOrbit = Orbites(a_i*1e3, e_i, 0, 0)./1e3;     % divide by 1e3 because orbites fucntion returns values in meters
TransferOrbit = Orbites(a_t*1e3, e_t, 0, 0)./1e3;


lon = rad2deg(theta_entry) - rad2deg(downrange);
lat = zeros(length(alt), 1);
lat = lat + 1;
alt = (alt./1e3 + Mars.radius);
Pos = SphereToCartesian(alt, rad2deg(theta_entry), lat);
T = SphereToCartesian(alt, lon, 90);
figure('color','white');
set(gcf, 'Position', [0 0 3840, 2160]);
grid on;
MarsSphere(gcf, 'km');
hold on
% plot3(InitialOrbit(361, 1), InitialOrbit(361, 2), InitialOrbit(361, 3), 'rx', 'LineWidth', 1.5)
plot3(InitialOrbit(180:361, 1), InitialOrbit(180:361, 2), InitialOrbit(180:361, 3), '-b', 'LineWidth', 1.5)
plot3(InitialOrbit(180, 1), InitialOrbit(180, 2), InitialOrbit(180, 3), 'go', 'Linewidth', 1.5)
plot3(TransferOrbit(60, 1), TransferOrbit(60, 2), TransferOrbit(60, 3), 'rx', 'Linewidth', 2)
plot3(TransferOrbit(60:181, 1), TransferOrbit(60:181, 2), TransferOrbit(60:181, 3), '-c', 'LineWidth', 1.5)
plot3(T(:, 1), T(:, 2), T(:, 3), '-y', 'LineWidth', 1.5)
ax = gca;
ax.Clipping = 'off';
hold off
xlabel('Pos x [km]');
ylabel('Pos y [km]');
zlabel('Pos z [km]');
grid on;
axis equal


burnTime = Vehicle.mass/Vehicle.MassFlow * (1-exp(-abs(DeltaV*1000)/Vehicle.Isp))   % s - duration of the deceleration burn
mass_afterBurn = Vehicle.mass - Vehicle.MassFlow * burnTime;                % kg - mass of Vehicle after burn

flight_path_angle = atan2(e_t * sin(theta_entry), 1 + e_t * cos(theta_entry));
fprintf('flight path angle = %.4f°\n', rad2deg(flight_path_angle));
V_entry = V_entry * 1e3;   % we need the value in m/s
save('Data/orbit.mat', 'V_entry', 'flight_path_angle');



