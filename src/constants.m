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

Vehicle.J = 928.951;
Vehicle.S = pi * ((Vehicle.Dimensions.d1 / 2) * cosd(Vehicle.Dimensions.theta1))^2 ; % m^2 - Vehicle aerodynamic surface
Vehicle.d = 0.53113125; % m - Vehicle aerodynamic dimension
Vehicle.CMalpha = -0.07; % Couple coefficient
Vehicle.CMq = -0.05; % Damping coefficent
Vehicle.CMdelta = 0.10; % High-lift device coefficient
Vehicle.Thrust = 400;   % N
Vehicle.Isp = 3000;     % m/s
Vehicle.MassFlow = Vehicle.Thrust / Vehicle.Isp; % kg/s

Vehicle.Parachute.CD = 1.75; % Drag coefficient
Vehicle.Parachute.D = 40; % m - Parachute diameter
Vehicle.Parachute.S =  4 * pi * Vehicle.Parachute.D^2 / 8; % m^2 - Surface area

%% Export constants
save('Data/constants.mat');
