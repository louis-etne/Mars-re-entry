clear all;
close all;
clc;

%% Mars constants
Mars.mass = 6.4185e23; % kg - The mass of Mars
Mars.radius = 3389.5; % km - Mean Mars radius
Mars.g0 =  3.711; % m/s^2 - Surface gravity
Mars.mu = 42828.37; % km^3/s^2 - Standard gravitational parameter

%% Mars atmosphere constants
Atm.hs = 11.1e3; % km - Atmospheric scale height
Atm.rho0 = 0.020; % kg/m^3 - Atmoshpere density at h = 0
Atm.gamma = 1.29;
Atm.R = 191.8;

%% Orbit constants
Orbit.height = 500; % km - Orbit height
Orbit.altitude = Mars.radius + Orbit.height; % km - Orbit altitude

%% Vehicle constants
Vehicle.mass = 50; % kg - Vehicle mass
Vehicle.J = 1.5; % kg/m^2 - Vehicle intertia
Vehicle.S = 0.80; % m^2 - Vehicle aerodynamic surface
Vehicle.d = 0.05; % m - Vehicle aerodynamic dimension
Vehicle.CD0 = 1.20; % Drag coefficient
Vehicle.CLalpha = 0.80; % Lift coefficient
Vehicle.CMalpha = -0.07; % Couple coefficient
Vehicle.CMq = -0.05; % Damping coefficient
Vehicle.CMdelta = 0.10; % Lift device coefficient

%% Export constants
save('Data/constants.mat');