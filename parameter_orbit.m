clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat', 'r', 'radius_m', 'mu_m');

%% Parameter of orbit around mars : 
% HYPOTHESIS : 
%   - circular orbit at 500 km altitude

% speed departure orbit
V_departure_orbit = sqrt(2 * ( -mu_m / (2 * r) + mu_m / r));