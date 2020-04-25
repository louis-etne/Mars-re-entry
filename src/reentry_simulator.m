clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat');

%% Initial conditions
v_ini = 3627.7; %  m/s - Initial velocity at h_ini
gamma_ini = deg2rad(-20.5); % rad - Angle between horizontal plane and the velocity vector - Initial flight path angle
h_ini = 60000; % m - Initial altitude
phi_ini = deg2rad(0.0); % rad - Downrange

xi = [v_ini, gamma_ini, h_ini, phi_ini];
ts = [0 125]; % s - Time span, [t0 tf]

%% Simulation
options = odeset('reltol', 1e-8);
[t, x0] = ode45(@(t0,y0)Dynamics(t0, y0, Mars, Atm, Vehicle, deg2rad(0)), ts, xi, options);
x0 = real(x0);

[tm10, xm10] = ode45(@(tm10,ym10)Dynamics(tm10, ym10, Mars, Atm, Vehicle, deg2rad(-10)), ts, xi, options);
xm10 = real(xm10);

[t10, x10] = ode45(@(t10,y10)Dynamics(t10, y10, Mars, Atm, Vehicle, deg2rad(10)), ts, xi, options);
x10 = real(x10);

[t20, x20] = ode45(@(t20,y20)Dynamics(t20, y20, Mars, Atm, Vehicle, deg2rad(20)), ts, xi, options);
x20 = real(x20);

%% Assignations
v0 = x0(:, 1);
h0 = x0(:, 3);

v10 = x10(:, 1);
h10 = x10(:, 3);

vm10 = xm10(:, 1);
hm10 = xm10(:, 3);

v20 = x20(:, 1);
h20 = x20(:, 3);

figure;
hold on;
plot(t, h0);
plot(tm10, hm10);
plot(t10, h10);
plot(t20, h20);
legend('\alpha = 0°', '\alpha = -10°', '\alpha = 10°', '\alpha = 20°');
grid on;
title('Altitude vs time');
xlabel('Time (s)');
ylabel('Altitude (m)');

figure;
hold on;
plot(t, v0);
plot(tm10, vm10);
plot(t10, v10);
plot(t20, v20);
legend('\alpha = 0°', '\alpha = -10°', '\alpha = 10°', '\alpha = 20°');
grid on;
title('Velocity vs time');
xlabel('Time (s)');
ylabel('Velocity (m.s^{-1})');

figure;
ax = axes;
hold on;
plot(h0, v0);
plot(hm10, vm10);
plot(h10, v10);
plot(h20, v20);
legend('\alpha = 0°', '\alpha = -10°', '\alpha = 10°', '\alpha = 20°');
grid on;
title('Velocity versus altitude');
xlabel('Altitude (m)');
ylabel('Velocity (m.s^{-1})');
set(ax, 'XDir', 'reverse');
