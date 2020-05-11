clear all; clc; close all;

addpath('Data');
addpath('Functions');

load('orbit.mat');
load('fit_temp.mat');
load('constants.mat');

aero_coefs.CA = load('fit_areo_axial_coef.mat');
aero_coefs.CN = load('fit_areo_norm_coef.mat');



%% Initial conditions
v_ini = V_entry; %  m/s - Initial velocity at h_ini
gamma_ini = -flight_path_angle; % rad - Angle between horizontal plane and the velocity vector - Initial flight path angle
h_ini = 120000; % m - Initial altitude
phi_ini = deg2rad(0.0); % rad - Downrange
theta_ini = deg2rad(-80); % rad - Angle between horizontal and the capsule axis
q_ini = deg2rad(0.0); % rad/s - Angular speed between horizontal and the capsule axis 

xi = [v_ini, gamma_ini, h_ini, phi_ini, theta_ini, q_ini];

%% Simulation parameters
sim_params.parachute_at = 3000; % m - Open parachute at 8km
sim_params.parachute_mach = 1.8; % Max mach number for parachute
sim_params.r_des = sim_params.parachute_at + Mars.radius * 1000;
sim_params.v_des = sim_params.parachute_mach * SpeedOfSound(sim_params.parachute_at, Atm);
sim_params.rho_des = Density(sim_params.parachute_at, Atm);

% Controller parameters
sim_params.tau = 2;
sim_params.zeta = 7;
sim_params.wn = 20;

sim_params.controlled = true;