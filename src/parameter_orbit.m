clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat', 'Mars', 'Orbit', 'Vehicle');


%% Parameter of orbit around mars : 
% HYPOTHESIS : 
%   - circular orbit at 500 km altitude
% Supposition for the calculation : Mars Atmosphere begin at 120 km high
% speed departure orbit
V_departure_orbit = sqrt(2 * (-Mars.mu / (2 * Orbit.radius) + Mars.mu / Orbit.radius));

r_p = Mars.radius + 60; % km - altitude of transfer orbit periapsis
a_i = Orbit.radius;      % km - semi-major axis of initial orbit
e_i = 0;                 % eccentricity of initial orbit
InitialOrbit = orbites(a_i*1e3, e_i, 0, 0)./1e3;     % divide by 1e3 because orbites fucntion returns values in meters

a_t = (Orbit.radius + r_p)/2;                        % km - sma of transfer orbit
e_t = (Orbit.radius - r_p) / (Orbit.radius + r_p);   % eccentricity of transfer orbit
Va_t = sqrt(Mars.mu * (2/Orbit.radius - 1/a_t));     % km/s velocity at apoapsis of transfer orbit
Vp_t = sqrt(Mars.mu * (2/r_p - 1/a_t));              % km/s velocity at periapsis of transfer orbit
DeltaV = V_departure_orbit - Va_t;

InitialOrbit = orbites(a_i*1e3, e_i, 0, 0)./1e3;         % divide by 1e3 because orbites fucntion returns values in meters
TransferOrbit = orbites(a_t*1e3, e_t, 0, 0)./1e3;

figure('color','white');
set(gcf, 'Position', [0 0 1920, 1080]);
grid on;
MarsSphere(gcf, 'km')
hold on
plot3(InitialOrbit(:, 1), InitialOrbit(:, 2), InitialOrbit(:, 3), '-b', 'LineWidth', 1.5)
plot3(TransferOrbit(:, 1), TransferOrbit(:, 2), TransferOrbit(:, 3), '-c', 'LineWidth', 1.5)
hold off
xlabel('Pos x [km]');
ylabel('Pos y [km]');
zlabel('Pos z [km]');
grid on;
axis equal

burnTime = Vehicle.mass/Vehicle.MassFlow * (1-exp(-DeltaV/Vehicle.Isp));    % s - duration of the deceleration burn
mass_afterBurn = Vehicle.mass - Vehicle.MassFlow * burnTime;                % kg - mass of Vehicle after burn




