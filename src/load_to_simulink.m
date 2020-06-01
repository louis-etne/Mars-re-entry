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
gamma_ini = 0.0587; % rad - Angle between horizontal plane and the velocity vector - Initial flight path angle
h_ini = 120000; % m - Initial altitude
phi_ini = deg2rad(0.0); % rad - Downrange
theta_ini = deg2rad(-80); % rad - Angle between horizontal and the capsule axis
q_ini = deg2rad(0.0); % rad/s - Angular speed between horizontal and the capsule axis 

xi = [v_ini, gamma_ini, h_ini, phi_ini, theta_ini, q_ini];

%% Simulation parameters
sim_params.parachute_at = 4000; % m - Open parachute at 8km
sim_params.parachute_mach = 1.8; % Max mach number for parachute
sim_params.r_des = sim_params.parachute_at + Mars.radius * 1000;
sim_params.v_des = sim_params.parachute_mach * SpeedOfSound(sim_params.parachute_at, Atm);
sim_params.rho_des = Density(sim_params.parachute_at, Atm);

% Controller parameters
sim_params.tau = 1.5;
sim_params.zeta = 10;
sim_params.wn = 20;

sim_params.controlled = true;


output = sim("prof_eq.slx", "RelTol", "1e-6", "StartTime", "0", "StopTime", "1000");


t = output.h.time;
h = output.h.data;
v = output.v.data;
a = output.a.data;
M = output.M.data;
gamma = output.gamma.data;
D = output.Daero.data;
L = output.Laero.data;

figure();
plot(t, h);
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude function of time");
figure();
plot(t, v);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title("Velocity function of time");
figure();
plot(h, v);
xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude");
set(gca,'Xdir','reverse');
figure();
plot(t, a);
title("Acceleration function of time");
xlabel("Time (s)");
ylabel("Acceleration");
figure();
plot(t, M);
title("Mach number function of time");
xlabel("Time (s)");
ylabel("Mach number");
figure();
plot(t, rad2deg(gamma));
title("Flight path angle \gamma function of time");
xlabel("Time (s)");
ylabel("Flight path angle (deg)");
figure();
hold on;
plot(t, D);
plot(t, L);
legend("Laero", "Daero");
xlabel("Time (s)");
ylabel("Aerodynamic forces");
title("L_{aero} and D_{aero} function of time");
