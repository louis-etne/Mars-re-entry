clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat');

%% Initial conditions
h_ini = 120000; % m - Initial altitude
v_ini = 6100; % m/s - Initial velocity at h_ini
s_ini = deg2rad(0); % rad - Downrange - Horizontal position
gamma_ini = deg2rad(-20.5); % rad - Angle between horizontal plane and the velocity vector - Initial flight path angle
theta_ini = deg2rad(-80); % rad - Angle between horizontal plane and the vehicle axis
q_ini = deg2rad(0); % rad - Pitch angular velocty : q = thetad

x0 = [h_ini, v_ini, s_ini, gamma_ini, theta_ini, q_ini];
ts = [0 125]; % s - Time span, [t0 tf]

%% Simulation
options = odeset('reltol', 1e-8);
[t, x] = ode45(@(t,y)Dynamics(t, y, Mars, Atm, Vehicle), ts, x0, options);
x = real(x);

%% Assignations
h = x(:, 1);
v = x(:, 2);
s = x(:, 3);
gamma = x(:, 4);
theta = x(:, 5);
q = x(:, 6);

figure;
plot(t, h);
title('Altitude vs time');
xlabel('Time (s)');
ylabel('Altitude (m)');

figure;
plot(t, v);
title('Velocity vs time');
xlabel('Time (s)');
ylabel('Velocity (m.s^{-1})');

figure;
ax = axes;
plot(h, v);
title('Velocity versus altitude');
xlabel('Altitude (m)');
ylabel('Velocity (m.s^{-1})');
set(ax, 'XDir', 'reverse');
