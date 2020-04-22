clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

%% Plot Mars atmosphere density
h = 0:120000; % m - Altitudes at which the density will be calculated
[T, P, rho1, rho2] = MarsAtmosphere(h);

%% Plots
figure;
hold on;
plot(rho, h);
ylabel('Altitude (m)');
xlabel('Density \rho (kg.m^{-3})');
title('Density versus Altitude');

figure;
hold on;
plot(T, h);
ylabel('Altitude (m)');
xlabel('Temperature (°C)');
title('Temperature versus Altitude');

figure;
hold on;
plot(P, h);
ylabel('Altitude (m)');
xlabel('Pressure (Pa)');
title('Pressure versus Altitude');

figure;
ax = axes;
hold on;
plot(T, rho);
xlabel('Temperature (°C)');
ylabel('Density \rho (kg.m^{-3})');
title('Density versus Temperature');
set(ax, 'Xdir', 'reverse');

figure;
ax = axes;
hold on;
plot(P, rho);
xlabel('Pressure (Pa)');
ylabel('Density \rho (kg.m^{-3})');
title('Density versus Pressure');
set(ax, 'Xdir', 'reverse');