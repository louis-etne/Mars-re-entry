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

r_p = Mars.radius + 0;   % km - altitude of transfer orbit periapsis
a_i = Orbit.radius;      % km - semi-major axis of initial orbit
e_i = 0;                 % eccentricity of initial orbit

a_t = (Orbit.radius + r_p)/2;                        % km - sma of transfer orbit
e_t = (Orbit.radius - r_p) / (Orbit.radius + r_p);   % eccentricity of transfer orbit
Va_t = sqrt(Mars.mu * (2/Orbit.radius - 1/a_t));     % km/s velocity at apoapsis of transfer orbit
Vp_t = sqrt(Mars.mu * (2/r_p - 1/a_t));              % km/s velocity at periapsis of transfer orbit
DeltaV = Va_t - V_departure_orbit;
h = r_p * Vp_t;
p = h^2 / Mars.mu;
r_entry = Mars.radius + 120;
theta_entry = acos(p/(r_entry * e_t) - 1/e_t);
V_entry = sqrt(Mars.mu * (2/r_entry - 1/a_t));

InitialOrbit = flip(Orbites(a_i*1e3, e_i, 0, 0)./1e3);     % divide by 1e3 because orbites fucntion returns values in meters
TransferOrbit = flip(Orbites(a_t*1e3, e_t, 0, 0)./1e3);

lon = rad2deg(theta_entry) - rad2deg(downrange);
lat = zeros(length(alt), 1);
lat = lat + 1;
alt = (alt./1e3 + Mars.radius);
Pos = SphereToCartesian(alt, lon, 90);
figure('color','white');
set(gcf, 'Position', [0 0 3840, 2160]);
grid on;
MarsSphere(gcf, 'km');
xlabel('Pos x [km]');
ylabel('Pos y [km]');
zlabel('Pos z [km]');
axlen = Mars.radius*1.2;
axis([-axlen, axlen, -axlen, axlen, -axlen, axlen])
ax = gca;
ax.Clipping = 'off';
index1 = int16(length(InitialOrbit)/2);
index2 = int16(length(TransferOrbit)/2);

hold on
%%
% ellipse1 = animatedline('Color', 'b', 'LineWidth', 1.5);
% ellipse2 = animatedline('Color', 'c', 'LineWidth', 1.5);
% descent = animatedline('Color', 'g', 'LineWidth', 2);

% for i = 1:int16(length(InitialOrbit)/2)
%    addpoints(ellipse1, InitialOrbit(i, 1), InitialOrbit(i, 2), InitialOrbit(i, 3));
%    drawnow;
% end
% for i = int16(length(TransferOrbit)/2):3020
%    addpoints(ellipse2, TransferOrbit(i, 1), TransferOrbit(i, 2), TransferOrbit(i, 3));
%    drawnow;
% end
% for i = 1:length(Pos)
%     addpoints(descent, Pos(i, 1), Pos(i, 2), Pos(i, 3));
%     drawnow;
% end
%%
% comet3(ax, InitialOrbit(1:index1, 1), InitialOrbit(1:index1, 2), InitialOrbit(1:index1, 3))%, '-b', 'LineWidth', 1.5)
% comet3(ax, TransferOrbit(index2:3020, 1), TransferOrbit(index2:3020, 2), TransferOrbit(index2:3020, 3))%, '-c', 'LineWidth', 1.5)
% comet3(ax, Pos(:, 1), Pos(:, 2), Pos(:, 3))%, '-y', 'LineWidth', 1.5)
%%

trajectory = [InitialOrbit(1:index1, :); TransferOrbit(index2:302, :); Pos(:, :)];
indexEnd = length(trajectory);
CustomComet3(ax, trajectory(:, 1), trajectory(:, 2), trajectory(:, 3))
% plot3(trajectory(indexEnd, 1), trajectory(indexEnd, 2), trajectory(indexEnd, 3), 'gx', 'markersize', 3);

burnTime = Vehicle.mass/Vehicle.MassFlow * (1-exp(-abs(DeltaV)/Vehicle.Isp));    % s - duration of the deceleration burn
mass_afterBurn = Vehicle.mass - Vehicle.MassFlow * burnTime;                % kg - mass of Vehicle after burn

flight_path_angle = atan2(e_t * sin(theta_entry), 1 + e_t * cos(theta_entry));
fprintf('flight path angle = %.4f°\n', rad2deg(flight_path_angle));
V_entry = V_entry * 1e3;   % we need the value in m/s
save('Data/orbit.mat', 'V_entry', 'flight_path_angle');



