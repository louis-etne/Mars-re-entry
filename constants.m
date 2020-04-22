clear all;
close all;
clc;

%% Mars constants
mars_mass = 6.4185e23; % kg - The mass of Mars
mars_volume = 163.18e9; % km^3 - The volume of Mars
mars_mean_density = 3933.5; % kg/m^3 - The mean density of Mars
g0 = 3.711; % m/s^2 - Surface gravity

%% Mars atmosphere constants
hs = 11.1e3; % km - Atmospheric scale height
rho0 = 0.020; % kg/m^3 - Atmoshpere density at h = 0

%% Export constants
save('Data/constants.mat');