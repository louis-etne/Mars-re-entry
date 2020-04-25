clear all;
close all;
clc;

%% Physical constants
Physics.R = 8.31446261815324; % J/(mol.K) - Universal Gas Constant

%% Mars constants
Mars.mass = 6.4185e23; % kg - The mass of Mars
Mars.radius = 3389.5; % km - Mean Mars radius
Mars.g0 =  3.711; % m/s^2 - Surface gravity
Mars.mu = 42828.37; % km^3/s^2 - Standard gravitational parameter

%% Mars atmosphere constants
Atm.hs = 11.1e3; % km - Atmospheric scale height
Atm.rho0 = 0.020; % kg/m^3 - Atmoshpere density at h = 0
Atm.gamma = 1.29; % Heat capacity ratio
Atm.mean_molar_mass = 43.34e-3; % kg/mol - Mean molar mass of gas
Atm.R = Physics.R / Atm.mean_molar_mass; % K/(kg.K) - Specific gas constant

%% Orbit constants
Orbit.altitude = 500; % km - Orbit altitude
Orbit.radius = Mars.radius + Orbit.altitude; % km - Orbit radius

%% Vehicle constants
Vehicle.mass = 1500; % kg - Vehicle mass

% Dimensions - See Fig. 4.91 of Aerodynamic Data of Space Vehicles (p. 114)
Vehicle.Dimensions.d1 = 3.5005; % m
Vehicle.Dimensions.d2 = 0.7623; % m
Vehicle.Dimensions.l1 = 1.6438; % m
Vehicle.Dimensions.l2 = 0.7500; % m
Vehicle.Dimensions.l3 = 0.8061; % m
Vehicle.Dimensions.r1 = 0.8762; % m
Vehicle.Dimensions.r2 = 0.0262; % m
Vehicle.Dimensions.theta1 = 20.00; % degrees
Vehicle.Dimensions.theta2 = 40.00; % degrees
Vehicle.Dimensions.theta3 = 62.18; % degrees

Vehicle.S = pi * ((Vehicle.Dimensions.d1 / 2) * cosd(Vehicle.Dimensions.theta1))^2 ; % m^2 - Vehicle aerodynamic surface
Vehicle.d = 0.05; % m - Vehicle aerodynamic dimension
Vehicle.C_A = 1.5; % Axial force coefficient
Vehicle.C_N = 0.03; % Normal force coefficient

%% Export constants
save('Data/constants.mat');
