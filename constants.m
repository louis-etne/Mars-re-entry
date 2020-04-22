clear all;
close all;
clc;

%% Mars constants
mass_m = 6.4185e23; % kg - The mass of Mars
volume_m = 163.18e9; % km^3 - The volume of Mars
mean_density_m = 3933.5; % kg/m^3 - The mean density of Mars
g0 = 3.711; % m/s^2 - Surface gravity
radius_m = 3389.5;
mu_m = 42828;

%% Mars atmosphere constants
hs = 11.1e3; % km - Atmospheric scale height
rho0 = 0.020; % kg/m^3 - Atmoshpere density at h = 0
gamma_m = 1.29;
R_m = 191.8;

%% Orbit constants
alt = 500; % km - Orbit altitude

% circular orbit
r = radius_m + alt;

%% Export constants
save('Data/constants.mat');