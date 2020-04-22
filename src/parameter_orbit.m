clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat', 'Mars', 'Orbit');

%% Parameter of orbit around mars : 
% HYPOTHESIS : 
%   - circular orbit at 500 km altitude

% speed departure orbit
V_departure_orbit = sqrt(2 * (-Mars.mu / (2 * Orbit.altitude) + Mars.mu / Orbit.altitude));